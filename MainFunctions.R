source('SubFunctions.R')

# ------------------------------------------------------------- #
# SPR_Imputation : Imputation step of SPR                       #
# ------------------------------------------------------------- #

## Input ------------------------------------------------------
# data: data with censored observations
# y_var: column name of observed survival times
# status_var: column name of censoring indicator
# K: the number of neighbors in KNN to impute LTF censored survival times (default = 10)
# EoR_distn: distribution to impute EoR censored survival times
# EoR_seed: seed number for sampling to impute EoR censored survival times (default = 1886)
## ------------------------------------------------------------

## Output -----------------------------------------------------
# $data: imputed data
# $LTF.idx: the index vector of LTF censored survival times
# $EoR.idx: the index vector of EoR censored survival times
# $EoR.distn.params: the parameters of distribution to impute EoR censored survival times
## ------------------------------------------------------------

k.rf = readRDS('modeling_k_rf.rds')   # for estimation of flattening parameter k
SPR_Imputation = function(data, y_var, status_var, K, EoR_distn, EoR_seed, EoR_time){
  data = data.frame(data)
  
  ## ----------------------------- Consider EoR_time ----------------------------- ##
  # death_index, censored_index, EoR_index
  death_index = which(data[,status_var]==1)
  censored_index = which(data[,status_var]==0)
  
  EoR.idx = which(data[, y_var] > EoR_time); 
  data[EoR.idx, status_var] = 0; data[EoR.idx, y_var] = EoR_time
  censored_index = which(data[,status_var]==0)
  LTF.idx = setdiff(censored_index, EoR.idx)
  
  ## -------------------- LTF censored Imputation and Flatten -------------------- ##
  # KNN Imputation
  Imp_dat = y_Imputation(data, y_var, status_var, n.cluster=K)
  EoR.idx = which(is.na(Imp_dat[,y_var]))
  # EoR.idx = which(data$status == 0 & data$eventtime >= EoR_time)
  Imp.idx = setdiff(which(data[,status_var] == 0), EoR.idx) #
  
  # flatten
  pred_dat = flat_k_data(NULL, data, Imp_dat, y_var, status_var, mean(data[,status_var] == 0))
  pred_dat = data.frame(t(pred_dat))
  k_to_flatten = predict(k.rf, pred_dat)$predictions
  flat_dat = Imp_flatten(Imp_dat, data, Imp.idx, EoR.idx, y_var, k = k_to_flatten)
  flat_dat[which(flat_dat$eventtime < data$eventtime),y_var] = data[which(flat_dat$eventtime < data$eventtime),y_var]
  
  ## -------------------------- EoR censored Imputation -------------------------- ##
  if(length(EoR.idx) != 0){
    RC_Imputation_D = flat_dat[setdiff(1:nrow(Imp_dat), EoR.idx),]
    EoR_D = Imp_dat[EoR.idx,]   # for Prediction
    
    # Linear Regression to estimate hazards ranking
    EoR_lm = lm(paste0(y_var, " ~.-", status_var), data = RC_Imputation_D)
    Reg_Fit = predict(EoR_lm, newdata = EoR_D)
    
    EoR_rank = rank(Reg_Fit)
    EoRmin = min(data[EoR.idx,y_var])
    
    # Imputation by Sampling
    ## Estimation of paramters
    if(toupper(EoR_distn) == 'NORMAL') {CLL = N_CLL} else
      if(toupper(EoR_distn) == 'WEIBULL') {CLL = Wei_CLL} else
        if(toupper(EoR_distn) == 'LOGNORMAL') {CLL = LN_CLL} else
        {return(print('Choose a distribution among Normal, Weibull and lognormal'))}
    if(toupper(EoR_distn) == 'NORMAL') {init = c(mean(data[,y_var]), sd(data[,y_var]))} else
      if(toupper(EoR_distn) == 'WEIBULL') {init = c(shape = 1, scale = 1)} else
        if(toupper(EoR_distn) == 'LOGNORMAL') {init = c(mean(data[,y_var]), sd(data[,y_var]))}
    # init = c(mean(data[,y_var]), sd(data[,y_var]))
    E_Dist = optim(init, function(x) -CLL(x, data[,y_var], data[,status_var]))
    
    ## Sampling
    sstrain = EoR_Imputation(data = flat_dat, y_var, status_var, y_rank = EoR_rank, EoR = EoRmin, dist = EoR_distn, E_Dist$par, EoR_seed = EoR_seed)
    out_list = list(data = sstrain, LTF.idx = Imp.idx, EoR.idx = EoR.idx, EoR.distn.params = E_Dist$par)
  }else{
    sstrain = flat_dat
    out_list = list(data = sstrain, LTF.idx = Imp.idx, EoR.idx = EoR.idx, EoR.distn.params = NULL)
  }
  
  return(out_list)
}


# ------------------------------------------------------------- #
# SPR_ASP : Calculation of ASP for SPR                          #
# ------------------------------------------------------------- #

## Input ------------------------------------------------------
# SPR_Surv: List of survival functions for newdata observations
# newdata: data to calculate ASP
# y_var: column name of observed survival times
# status_var: column name of censoring indicator
# max_time: Maximum survival time of survival function
## ------------------------------------------------------------

## Output -----------------------------------------------------
# ASP vector (for 1, 0.9,...,0.1, 0)
## ------------------------------------------------------------

SPR_ASP = function(SPR_Surv, newdata, y_var, status_var, max_time, EoR_time){
  newdata = data.frame(newdata)
  # Find survival times with survival probability of 0, 0.1, ..., 0.9, 1
  Surv_time = as.data.frame(t(sapply(1:nrow(newdata), function(x) Check_My_Surv(SPR_Surv[[x]], max_time, seq(1, 0.1, -0.1), EoR_time))))
  colnames(Surv_time) = paste0('P', seq(1, 0.1, -0.1))
  # Calculate ASP
  d = cal_ASP(Surv_time, newdata[,y_var], newdata[,status_var])
  return(d)
}

# ------------------------------------------------------------- #
# Tune_sd : Tuning sd for calculating ASP of SPR                #
# ------------------------------------------------------------- #

## Input ------------------------------------------------------
# sd_cand: Candidates of standard deviation for sampling continuous variables (default = from 0.5 to 2 by 0.1 sequence)
# fitted_model: fitted regression model using imputed data
# newdata: data to calculate ASP
# y_var: column name of observed survival times
# status_var: column name of censoring indicator
# conti_idx: column index vector of continuous variables on the newdata
# discr_idx: column index vector of discrete variables on the newdata
# num_sample: the number of samples to estimate survival function for each observation
# sample_seed: seed number for sampling (default = 54)
# max_time : Maximum survival time of survival function
# min_time: minimum survival time of survival function (default = NULL)
# BCRF_coef: Coefficients for BCRF when regression method is BCRF (default = NULL)
## ------------------------------------------------------------

## Output -----------------------------------------------------
# $SurvF: estimated survival function using the optimal sd for each observation (newdata)
# $ASP: optimal ASP vector
# $opt_sd: optimal standard deviation for sampling continuous variables
## ------------------------------------------------------------

Tune_sd = function(sd_cand, fitted_model, newdata, y_var, status_var, conti_idx, discr_idx,
                   num_sample, sample_seed=54,
                   max_time, min_time=NULL, EoR_time, BCRF_coef=NULL){
  
  if(class(fitted_model) == "ranger") {Surv_ftn = RF_Surv_ftn} else
    if(class(fitted_model) == "lm") {Surv_ftn = Lin_Surv_ftn} else
      if(class(fitted_model) == "xgb.Booster") {Surv_ftn = XGB_Surv_ftn}
  
  # Find the optimal standard deviation for sampling continuous variables
  ASP_res = list()
  for(cc in 1:length(sd_cand)){
    
    SPR_Surv = Surv_ftn(conti_range=sd_cand[cc], fitted_model=fitted_model, newdata=newdata,y_var, conti_idx=conti_idx, discr_idx=discr_idx, num_sample = num_sample, sample_seed = sample_seed, max_time=max_time, min_time=min_time, BCRF_coef=BCRF_coef)
    Surv_time = as.data.frame(t(sapply(1:nrow(newdata), function(x) Check_My_Surv(SPR_Surv[[x]], max_time, seq(1,0.1,-0.1), EoR_time))))
    colnames(Surv_time) = paste0('P', seq(1, 0.1, -0.1))
    rownames(Surv_time) = paste0("Obs",1:nrow(newdata))

    ASP_ = SPR_ASP(SPR_Surv, newdata, y_var, status_var, max_time, EoR_time)
    ASP_res[[cc]] = list(SurvF = SPR_Surv, ASP = ASP_, Surv_time = Surv_time)
    cat_asp = sum((ASP_res[[cc]]$ASP - seq(1, 0.1, -0.1))^2)
    cat('sd_cand: ', round(sd_cand[cc],1),'\n')
    # cat('sd_cand: ', round(sd_cand[cc],1) , ', ASP:', round(cat_asp,4),'\n')
  }

  opt_sd = which.min(lapply(lapply(ASP_res, '[[', 2), function(asp) mean((seq(1, 0.1, -0.1)-asp)^2)))
  opt_ASP = ASP_res[[opt_sd]]; opt_ASP$opt_sd = sd_cand[opt_sd]; opt_ASP$opt_asp = ASP_res[[opt_sd]]$ASP
  return(opt_ASP)
}


# ------------------------------------------------------------- #
# Cox_ASP : Calculation of ASP for Cox PH                       #
# ------------------------------------------------------------- #

## Input ------------------------------------------------------
# Cox_model: fitted CoxPH model
# newdata: data to calculate ASP
# y_var: column name of observed survival times
# status_var: column name of censoring indicator
# max_time: Maximum survival time of survival function
## ------------------------------------------------------------

## Output -----------------------------------------------------
# ASP vector (for 1, 0.9,...,0.1, 0)
## ------------------------------------------------------------

Cox_ASP = function(Cox_model, newdata, y_var, status_var, max_time, EoR_time){
  newdata = data.frame(newdata)
  # Estimation of survival function for each observation
  a = lapply(1:nrow(newdata), function(x) survfit(Cox_model, newdata = newdata[x,]))
  b = lapply(1:nrow(newdata), function(x) SurvF(data.frame(time = a[[x]]$time, surv = a[[x]]$surv), max_time))
  # Find survival times with survival probability of 0, 0.1, ..., 0.9, 1
  c = as.data.frame(t(sapply(1:nrow(newdata), function(x) Check_My_Surv(b[[x]], max_time, seq(1, 0.1, -0.1), EoR_time))))
  colnames(c) = paste0('P', seq(1, 0.1, -0.1))
  # Calculate ASP
  d = cal_ASP(c, newdata[,y_var], newdata[,status_var])
  return(d)
}

# ------------------------------------------------------------- #
# RSF_ASP : Calculation of ASP for RSF                          #
# ------------------------------------------------------------- #

## Input ------------------------------------------------------
# RSF_model: fitted RSF model
# newdata: data to calculate ASP
# y_var: column name of observed survival times
# status_var: column name of censoring indicator
# max_time: Maximum survival time of survival function
## ------------------------------------------------------------

## Output -----------------------------------------------------
# ASP vector (for 1, 0.9,...,0.1, 0)
## ------------------------------------------------------------

RSF_ASP = function(RSF_model, newdata, y_var, status_var, max_time, EoR_time){
  newdata = data.frame(newdata)
  # Estimation of survival function for each observation
  a = predict(RSF_model, newdata = newdata)
  b = lapply(1:nrow(newdata), function(x) SurvF(data.frame(time = a$time.interest, surv = a$survival[x, ]), max_time))
  # Find survival times with survival probability of 0, 0.1, ..., 0.9, 1
  c = as.data.frame(t(sapply(1:nrow(newdata), function(x) Check_My_Surv(b[[x]], max_time, seq(1, 0.1, -0.1), EoR_time))))
  colnames(c) = paste0('P', seq(1, 0.1, -0.1))
  # Calculate ASP
  d = cal_ASP(c, newdata[,y_var], newdata[,status_var])
  return(d)
}
