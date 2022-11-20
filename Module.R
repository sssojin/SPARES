# -----------------------------------------------------------------
# SPARES_Imputation

# Input
# - X : covariates matrix
# - y : survival times vector
# - censor_info : censoring vector
# - factor_var: vector of factor variables among independent variables
# - dist : distribution for imputation of censored survival times (weibull or lognormal)
# - reg_model : regression model for prediction of survival times (linear regression or random forests)
# - EoR_time : time of end of research
# - K : k value in KNN

# Output
# - y_hat : predicted survival times
# - imputed_y : imputed survival times for censored observation
# - model : fitted regression model
# - dat : new data with imputed_y
# -----------------------------------------------------------------

# SPARES_Imputation function
SPARES_Imputation = function(X, y, censor_info, dist, reg_model, EoR_time, K, 
                             y_var, status_var, discr_names, ...){
  ## -------------------------- 1. Setting dataset -------------------------- ##
  dat = data.frame(eventtime = y, status = censor_info, X);
  dat = dat %>% mutate_at(discr_names, as.factor) %>% as.data.frame()
  
  
  ## -------------------------- 2. Censored data Imputation -------------------------- ##
  EoR_seed_check = tryCatch({is.null(EoR_seed)}, error = function(e) print("EoR seed is randomly generated"))
  if(EoR_seed_check == "EoR seed is randomly generated") EoR_seed = as.integer(runif(1,0,999999))
  if(toupper(dist) == "WEIBULL"){
    sstrain = SPR_Imputation(data=dat, y_var=y_var, status_var=status_var, K=K, EoR_distn='weibull', EoR_seed=EoR_seed,  EoR_time = EoR_time)
  }else if(toupper(dist) == "LOGNORMAL"){
    sstrain = SPR_Imputation(data=dat, y_var=y_var, status_var=status_var, K=K, EoR_distn='lognormal', EoR_seed=EoR_seed, EoR_time = EoR_time)
  }else{
    print("dist is weibull or lognormal")
  }
  imputed_y = sstrain$data[,y_var]
  
  ## ---------------------------- 3. Fit Regression model ---------------------------- ##
  if(toupper(reg_model) == "LIN"){
    mdl = lm(eventtime ~., data = sstrain$data)
    y_hat = calib_survtime(pred_survtime=predict(mdl), min_survtime=min(dat$eventtime))
  }else if(toupper(reg_model) == "RF"){
    mdl = ranger(eventtime ~., data =  sstrain$data, num.trees=100, seed=EoR_seed)
    # y_hat = calib_survtime(pred_survtime=mdl$predictions, min_survtime=min(dat$eventtime))
    mdl_BCRF = predict.BC.SLR(mdl, sstrain$data, sstrain$data, 'eventtime'); BCRF_coef = mdl_BCRF$coef
    y_hat = calib_survtime(pred_survtime=mdl_BCRF$pred[,2], min_survtime=min(dat$eventtime))
  }else{
    print("reg_model is lin (linear regerssion) or rf (random forests)")
  }

  return(list(model = mdl, y_hat = y_hat, dat = dat, imputed_y = imputed_y))
}  
# -----------------------------------------------------------------
# SPARES_ASP

# Input
# - dat : output from the SPARES_Imputation
# - EoR_time : time of end of research

# Output
# - surv_time : survival probability matrix for each observation
# - asp : asp for the model
# -----------------------------------------------------------------

# SPARES_ASP function
SPARES_ASP = function(dat, y_var, status_var, discr_names,  EoR_time, ...){
  
  min_time = min(dat[,y_var])
  max_time = max(dat[,y_var])+sd(dat[,y_var])
  conti_idx = which(!(colnames(dat) %in% c(discr_names, y_var, status_var)))
  discr_idx = which(colnames(dat) %in% discr_names)
  
  survival_fn = Tune_sd(sd_cand=seq(1.5, 2.5, 0.1), fitted_model=model, newdata=dat, y_var=y_var, status_var=status_var,
                        conti_idx=conti_idx, discr_idx=discr_idx,
                        num_sample=500, max_time=max_time, min_time=min_time,
                        EoR_time = EoR_time, BCRF_coef=NULL)

  Opt_sd = survival_fn$opt_sd
  opt_ASP = survival_fn$opt_asp
  ASP = sum((opt_ASP-seq(1,0.1,-0.1))^2)
  Surv_time = survival_fn$Surv_time
  return(list(Surv_time = Surv_time, ASP = round(ASP,4), Opt_sd=Opt_sd, opt_ASP = opt_ASP))
}

