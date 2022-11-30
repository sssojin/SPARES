# Load pacakges
library(survival); library(randomForestSRC)
if(!require(ranger)) print("Installation of the package ranger recommended") else library(ranger)
if(!require(xgboost)) print("Installation of the package xgboost recommended") else library(xgboost)
library(dplyr); library(e1071); 
library(parallel); library(doParallel)
library(GeneralizedHyperbolic)
library(dplyr)
# library(GeneralizedHyperbolic)


# ------------------------------------------------------------- #
# LTF-Imputation ) y_Imputation, est_y                          #
# Flattening ) flat_k_data, Imp_flatten                         #
# EoR-Imputation ) N_CLL, Wei_CLL, LN_CLL, EoR_Imputation       #
# ------------------------------------------------------------- #
# Related Main Function : SPR_Imputation

# y_Imputation: KNN Imputation for LTF censored survival times
y_Imputation = function(data, y_var, status_var, n.cluster){
  rownames(data) = 1:nrow(data)
  censored_index = which(data[,status_var] == 0); death_index = which(data[,status_var] == 1); 
  y = data[,y_var]
  
  # KNN Imputation
  dissimilar = as.matrix(cluster::daisy(data %>% dplyr::select(-all_of(y_var), -status_var), metric = 'gower', stand = TRUE))
  estimates = sapply(censored_index, function(x) est_y(y, x, death_index, n.cluster, dissimilar))
  
  # Random Censoring 대체
  data[censored_index, y_var] = estimates
  return(data)
}

## est_y: Estimation of survival times by KNN (in y_Imputation)
est_y = function(y, myindex, d.index, n.cluster, dissimilar){
  # Complete observations who have longer survival times 
  larger_y = which(y >= y[myindex]); Lcandidate = intersect(larger_y, d.index)
  # Check if this censored one should be EoR Imputed
  if(length(Lcandidate) < n.cluster) return(NA)
  
  # KNN Imputation
  y.idx = dissimilar[myindex, d.index]
  pts.idx = as.numeric(names(sort(y.idx)[1:n.cluster]))
  imputed_cand = mean(y[pts.idx])
  
  if(imputed_cand < y[myindex]) return(y[myindex]) else return(imputed_cand)
}


# flat_k_data: To construct data for estimation of flattening parameter k
flat_k_data = function(Org_dat, Censor_dat, Imp_dat, y_var, status_var, censor){
  # Construct survival times for comparison
  EoR.idx = which(is.na(Imp_dat[,y_var])); Imp.idx = setdiff(which(Censor_dat[,status_var] == 0), EoR.idx)
  sim.Imp.y = sapply(1:nrow(Imp_dat), function(idx){
    if(idx %in% EoR.idx) return(Censor_dat[idx,y_var]) else return(Imp_dat[idx,y_var])
  })
  Imp_T = sim.Imp.y[Imp.idx]; Org_T = Org_dat[Imp.idx,y_var]
  
  # If Org_dat is not null
  if(is.null(Org_dat)) opt_k = NULL else opt_k = sd(Org_T)/sd(Imp_T)
  
  y_stats = apply(cbind(Censor_dat[,y_var], sim.Imp.y), 2, function(y){
    c(sd(y), diff(quantile(y, probs = c(0.25, 0.75))), e1071::skewness(y), e1071::kurtosis(y))
  })
  colnames(y_stats) = rownames(y_stats) = NULL
  return(c(k = opt_k,
           sd_ratio = y_stats[1,1]/y_stats[1,2], IQR_ratio = y_stats[2,1]/y_stats[2,2],
           skew_diff = diff(y_stats[3,]), kurt_diff = diff(y_stats[4,]),
           censor = censor))
}


# Imp_flatten : To flatten the KNN imputed survival times
Imp_flatten = function(Imp_dat, Censor_dat, Imp.idx, EoR.idx, y_var, k = NULL){
  min.time = min(Censor_dat[,y_var])
  sim.Imp.y = sapply(1:nrow(Imp_dat), function(idx){
    if(idx %in% EoR.idx) return(Censor_dat[idx,y_var]) else return(Imp_dat[idx,y_var])
  })
  Imp_T = sim.Imp.y[Imp.idx]
  
  if(is.null(k)) return('Enter flattening parameter k')
  flat_T = sapply(k*(Imp_T - mean(Imp_T)) + mean(Imp_T), function(x) max(x, min.time))
  
  Imp_dat[Imp.idx,y_var] = flat_T
  return(Imp_dat)
}


# Censored Loglikelihood Functions
## Wei_CLL : Weibull Censored Log-Likelihood
Wei_CLL = function(par, y, status){
  lik = status*dweibull(y, par[1], par[2], log = T) + (1-status)*pweibull(y, par[1], par[2], lower.tail = F, log.p = T)
  return(sum(lik))
}
## N_CLL : Normal Censored Log-Likelihood
N_CLL = function(par, y, status){
  lik = status*dnorm(y, par[1], par[2], log = T) + (1-status)*pnorm(y, par[1], par[2], lower.tail = F, log.p = T)
  return(sum(lik))
}
## LN_CLL : Lognormal Censored Log-Likelihood
LN_CLL = function(par, y, status){
  lik = status*dlnorm(y, par[1], par[2], log = T) + (1-status)*plnorm(y, par[1], par[2], lower.tail = F, log.p = T)
  return(sum(lik))
}


# EoR_Imputation : To impute EoR censored survival times by sampling from a parametric distribution
EoR_Imputation = function(data, y_var, status_var, y_rank, EoR, dist, par, 
                          EoR_seed, Nsample = 10^8){
  set.seed(EoR_seed)
  y = data[,y_var]
  if (toupper(dist) == 'NORMAL'){
    est = (dplyr::filter(data.frame(x = rnorm(Nsample, par[1], par[2])), x > EoR))[1:length(y_rank),]
  } else if (toupper(dist) == 'WEIBULL'){
    est = (dplyr::filter(data.frame(x = rweibull(Nsample, par[1], par[2])), x > EoR))[1:length(y_rank),]
  } else {
    est = (dplyr::filter(data.frame(x = rlnorm(Nsample, par[1], par[2])), x > EoR))[1:length(y_rank),]
  }
  na.idx = which(is.na(data[,y_var]==TRUE))
  data[na.idx, y_var] = sort(est)[y_rank]
  # data[is.na(y), y_var] = sort(est)[y_rank]
  return(data)
}


# ------------------------------------------------------------- #
# SurvF, Check_My_Surv, cal_ASP                                 #
# ------------------------------------------------------------- #
# Related Main Function : SPR_ASP, Cox_ASP, RSF_ASP

# SurvF : To approximate survival function by interpolation
SurvF = function(data, max.time){
  min.y = min(data$time); max.y = max(data$time)
  small.interval = ifelse(min.y > 2, 10^4, 10^2)
  big.interval = ifelse((max.time - max.y) > 2, 10^4, 10^2)
  small.y = data.frame(time = seq(0, min.y, length.out = small.interval), surv = 1)
  if(max.y == max.time) big.y = data.frame(time = max.time, surv = 0) else
    big.y = data.frame(time = seq(max.y, max.time, length.out = big.interval), surv = 0)
  mydata = rbind(small.y[-nrow(small.y),], data, big.y[-1,])
  return(approxfun(mydata$time, mydata$surv, rule = 1:2))
}


# Check_My_Surv : To find a root with a specific survival probability
Check_My_Surv = function(surv_ftn, max.root, check_prob, EoR_time){
  myftn = function(prob){
    get_root = function(x) surv_ftn(x) - prob
    root = uniroot(get_root, c(0, max.root))$root
    # EoR check
    if(!is.null(EoR_time) & root>EoR_time) return(EoR_time) else return(root)
  }
  time.pts = sapply(check_prob, myftn)
  return(time.pts)
}

cal_ASP = function(est.survtimes, obs.survtimes, status){
  P_SSR = apply(est.survtimes, 2, function(x){
    alive = as.numeric(x <= obs.survtimes)
    Ind = 1-as.numeric(x > obs.survtimes)*as.numeric(status == 0)
    return(alive*Ind)
  })
  risk_Ind = apply(est.survtimes, 2, function(x){
    as.numeric(x > obs.survtimes)*as.numeric(status == 0)
  })
  return(colSums(P_SSR)/(length(obs.survtimes)-colSums(risk_Ind)))
}


# ------------------------------------------------------------- #
# Make_Sample_Data                                              #
# Lin_Surv_ftn, RF_Surv_ftn, XGB_Surv_ftn                       #
# ------------------------------------------------------------- #
# Related Main Function : Tune_sd

# Make_Sample_Data : To generate similar samples with the observation
Make_Sample_Data = function(myind, newdata, conti_idx, discr_idx, num_sample, conti_range, sample_seed){
  alpha = 1
  set.seed(sample_seed)
  myobs = newdata[myind,]
  
  # Sampling continuous variables
  conti_sd = apply(newdata[,conti_idx], 2, function(x) sd(x, na.rm = T))
  conti_min = myobs[,conti_idx] - conti_range*conti_sd
  conti_max = myobs[,conti_idx] + conti_range*conti_sd
  conti_sample = as.data.frame(sapply(1:ncol(conti_min), function(x) runif(n = num_sample, conti_min[,x], conti_max[,x])))
  #conti_sample = as.data.frame(sapply(1:ncol(conti_min), function(x) rnorm(n = num_sample, myobs[,conti_idx][x], sd = conti_sd[x]*alpha))) # normal dist
  
  # Sampling discrete variables
  discr_sample = as.data.frame(lapply(newdata[myind,discr_idx], rep, num_sample))
  
  if(length(discr_idx) != 0) xvar = cbind(conti_sample, discr_sample) else xvar = conti_sample
  colnames(xvar) = c(colnames(newdata[conti_idx]), colnames(newdata[discr_idx]))
  return(xvar)
}


# Estimation of survival function in SPR algorithm
## By generating residuals
Est_Surv = function(obsd_y, pred_y, test_y, method,
                    num_sample, sample_seed,
                    max_time, min_time, d, ...){
  non_cen_ind = which(Train$status == 1)
  non_cen_obsd_y = obsd_y[non_cen_ind]
  non_cen_pred_y = pred_y[non_cen_ind]
  # assign function
  # res = obsd_y-pred_y # original ver.
  res = non_cen_obsd_y - non_cen_pred_y # new ver. use only non-censored data
  m = num_sample*10
  
  if(toupper(method) == "NORMAL"){
    loglik_fun = function(x) {-sum(dnorm(res, 0, x, log = T))}
    sample_fun = function(x) {rnorm(num_sample, 0, x)}
  }else if(toupper(method) == "T-DISTN"){
    loglik_fun = function(x) {-sum(dt(res, df = x, log = T))}
    # sample_fun = function(x) {rt(num_sample, df = x)}
    sample_fun = function(x) {rt(m, df = x)}
  }else if(toupper(method) == "NIG"){
    loglik_fun = function(x) {-sum(log(dnig(res, param = x)))}
    sample_fun = function(x) {rnig(num_sample, mu = x[1], delta = x[2], alpha = x[3])} #
  }else if(toupper(method) == "GAMMA"){
    res_t = abs(res-max(res)) # reversing residual
    loglik_fun = function(x) {-sum(dgamma(res_t, shape = x[1], scale = x[2]))}
    sample_fun = function(x) {rgamma(num_sample, shape = x[1], scale = x[2])}
  }else if(toupper(method) == "SN"){
    loglik_fun = function(x) {-sum(dsn(res, xi = x[1], omega = x[2], alpha = x[3], log = T))}
    sample_fun = function(x) {rsn(num_sample,  xi = x[1], omega = x[2], alpha = x[3])}
  }else{
    stop("method is one of normal, t-distn and NIG")
  }
  
  # Estimation of paramter(s)
  if(toupper(method) %in% c("NORMAL", "T-DISTN")){
    est_par = optimize(function(x) loglik_fun(x), c(10^(-8), 10^8))$minimum
  }else if(toupper(method == "NIG")){
    #est_par = optim(c(0,1,1,0), function(x) loglik_fun(x))$par ## sojin
    est_par = nigFit(res)$param
  }else if(toupper(method) == "GAMMA"){
    est_par = optim(c(2,2), function(x) loglik_fun(x))$par
  }else if(toupper(method) == "SN"){
    est_par = optim(c(0, 0.1, -4), function(x) loglik_fun(x))$par
  }
  
  # Estimation of survival function for each observation
  SurvFList = lapply(1:length(test_y), function(kk){
    # kk = 1
    set.seed(sample_seed+kk)
    #sampled_res = sample_fun(est_par)
    if(toupper(method) %in% c("NORMAL", "T-DISTN","NIG","SN")) {
      sampled_res = sample_fun(est_par)
      #sampled_res = sample_fun(1)
    }else if(toupper(method) == "GAMMA"){
      sampled_res = abs(sample_fun(est_par) - min(sample_fun(est_par)))
    }
    
    # hist(sampled_res)
    
    lm1.pop = test_y[kk] + sampled_res
    # hist(lm1.pop)
    
    mysam1 = function(x, pop, n, d){
      # mydist1 = 1/abs((x - pop)^3)
      # mydist1 = 1/(x - pop)^2
      mydist1 = 1/abs(x - pop)^d
      res = sample(pop, n, replace = T, prob = mydist1)
      return (res)
    }
    #kk=500
    predicted = mysam1(sort(test_y)[kk], lm1.pop, 1000, d)
    #hist(predicted)
    #plot(sort(predicted))
    predicted_ = ifelse(predicted < 0, min_time, predicted)
    # hist(predicted_)
    
    # dup_pred = rep(test_y[kk], num_sample)
    # predicted_ = ifelse(dup_pred + sampled_res < 0, min_time, dup_pred+sampled_res)
    
    # Predict survival function by Kaplan-Meier method
    KM_ = survfit(Surv(predicted_, rep(1, length(predicted_))) ~ 1)
    KM_surv = data.frame(time = KM_$time, surv = KM_$surv)
    
    # ggplot(aes(x = time, y = surv), data = KM_surv) + geom_line()
    
    return(SurvF(KM_surv, max_time))
    # return(KM_surv)
  })
  return(list(Surv = SurvFList, par = est_par))
}



## Lin_Surv_ftn : SPR-lin
#Lin_Surv_ftn = function(conti_range, fitted_model, newdata, y_var,
#                        num_sample, sample_seed,
#                        max_time, min_time, ...){
#  newdata = data.frame(newdata)
#  pred_y = predict(fitted_model, newdata = newdata)

# Estimation of sigma
#  res = newdata[,y_var]-pred_y
#  sig_hat = optimize(function(x) -sum(dnorm(res, 0, x, log = T)), c(10^(-8), 10^8))$minimum

# Estimation of survival function for each observation
#  SurvFList = lapply(1:nrow(newdata), function(i){
#    set.seed(sample_seed+i)
#    sampled_res = rnorm(num_sample, 0, sig_hat)
#    predicted_ = ifelse(pred_y+sampled_res < 0, min_time, pred_y+sampled_res)

# Predict survival function by Kaplan-Meier method
#    KM_ = survfit(Surv(predicted_, rep(1, length(predicted_))) ~ 1)
#    KM_surv = data.frame(time = KM_$time, surv = KM_$surv)
#    return(SurvF(KM_surv, max_time))
#  })
#  return(SurvFList)
#}
Lin_Surv_ftn = function(conti_range, fitted_model, newdata, yvar,conti_idx, discr_idx,
                        num_sample, sample_seed,
                        max_time, min_time, BCRF_coef, ...){
  SurvFList = lapply(1:nrow(newdata), function(i){
    sampled_data = Make_Sample_Data(i, newdata, conti_idx, discr_idx, num_sample, conti_range, sample_seed)
    if(!is.null(BCRF_coef)){
      # When BCRF is applied
      pred_surv_times = predict(fitted_model, newdata = sampled_data) #$predictions
      pred_surv_times = BCRF_coef[1]+BCRF_coef[2]*pred_surv_times
      pred_surv_times = ifelse(pred_surv_times < 0, min_time, pred_surv_times)
    }else{
      pred_surv_times = predict(fitted_model, newdata = sampled_data) #$predictions
    }
    
    # Predict survival function by Kaplan-Meier method
    KM_ = survfit(Surv(pred_surv_times, rep(1, num_sample)) ~ 1)
    KM_surv = data.frame(time = KM_$time, surv = KM_$surv)
    return(SurvF(KM_surv, max_time))
  })
  return(SurvFList)
}


## RF_Surv_ftn : SPR-rf(or bcrf)
RF_Surv_ftn = function(conti_range, fitted_model, newdata, yvar,conti_idx, discr_idx,
                       num_sample, sample_seed,
                       max_time, min_time, BCRF_coef, ...){
  SurvFList = lapply(1:nrow(newdata), function(i){
    sampled_data = Make_Sample_Data(i, newdata, conti_idx, discr_idx, num_sample, conti_range, sample_seed)
    if(!is.null(BCRF_coef)){
      # When BCRF is applied
      pred_surv_times = predict(fitted_model, data = sampled_data)$predictions
      pred_surv_times = BCRF_coef[1]+BCRF_coef[2]*pred_surv_times
      pred_surv_times = ifelse(pred_surv_times < 0, min_time, pred_surv_times)[1:num_sample]
    }else{
      pred_surv_times = predict(fitted_model, data = sampled_data)[1:num_sample]
    }
    
    # Predict survival function by Kaplan-Meier method
    KM_ = survfit(Surv(pred_surv_times, rep(1, num_sample)) ~ 1)
    KM_surv = data.frame(time = KM_$time, surv = KM_$surv)
    return(SurvF(KM_surv, max_time))
  })
  return(SurvFList)
}

## XGB_Surv_ftn : SPR-xgb
XGB_Surv_ftn = function(conti_range, fitted_model, newdata, conti_idx, discr_idx,
                        num_sample, sampled_seed,
                        max_time, min_time, ...){
  SurvFList = lapply(1:nrow(newdata), function(i) {
    sampled_data = Make_Boot_Data(i, newdata, conti_idx, discr_idx, num_sample, conti_range, sample_seed)
    sampled_data = model.matrix(~.-1, sampled_data); sampled_data = sampled_data[,fitted_model$feature_names]
    pred_surv_times = predict(fitted_model, newdata = sampled_data)
    pred_surv_times = ifelse(pred_surv_times < 0, min_time, pred_surv_times)
    
    # Predict survival function by Kaplan-Meier method
    KM_ = survfit(Surv(pred_surv_times, rep(1, num_sample)) ~ 1)
    KM_surv = data.frame(time = KM_$time, surv = KM_$surv)
    return(SurvF(KM_surv, max_time))
  })
  return(SurvFList)
}


# ------------------------------------------------------------- #
# calib_survtime
# -> to calibrate predicted survival times by SPR
# ------------------------------------------------------------- #

calib_survtime = function(pred_survtime, min_survtime){
  ifelse(pred_survtime < 0, min_survtime, pred_survtime)
}
