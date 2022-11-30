
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
  Imputed_dat = sstrain$data
  imputed_y = sstrain$data[,y_var]
  mdl_dat =  sstrain$data%>% select(-status_var)
  
  ## ---------------------------- 3. Fit Regression model ---------------------------- ##
  if(toupper(reg_model) == "LIN"){
    mdl = lm(eventtime ~., data = mdl_dat)
    y_hat = calib_survtime(pred_survtime=predict(mdl), min_survtime=min(dat$eventtime)); BCRF_coef = NULL
  }else if(toupper(reg_model) == "RF"){
    mdl = ranger(eventtime ~., data =  mdl_dat, num.trees=100, seed=EoR_seed)
    # y_hat = calib_survtime(pred_survtime=mdl$predictions, min_survtime=min(dat$eventtime))
    mdl_BCRF = predict.BC.SLR(mdl, mdl_dat, mdl_dat, 'eventtime'); BCRF_coef = mdl_BCRF$coef
    y_hat = calib_survtime(pred_survtime=mdl_BCRF$pred[,2], min_survtime=min(dat$eventtime))
  }else{
    print("reg_model is lin (linear regerssion) or rf (random forests)")
  }

  return(list(model = mdl, y_hat = y_hat, dat = dat, Imputed_dat=Imputed_dat, imputed_y = imputed_y, BCRF_coef = BCRF_coef))
}  

# SPARES_ASP function
SPARES_ASP = function(dat, Imputed_dat, y_var, status_var, discr_names,  EoR_time, ...){
  
  min_time = min(dat[,y_var])
  max_time = max(dat[,y_var])+sd(dat[,y_var])
  conti_idx = which(!(colnames(Imputed_dat) %in% c(discr_names, y_var, status_var)))
  discr_idx = which(colnames(dat) %in% discr_names)
  
  survival_fn = Tune_sd(sd_cand=seq(1.5, 2.5, 0.1), fitted_model=model, newdata=Imputed_dat, y_var=y_var, status_var=status_var,
                        conti_idx=conti_idx, discr_idx=discr_idx,
                        num_sample=500, max_time=max_time, min_time=min_time,
                        EoR_time = EoR_time, BCRF_coef=BCRF_coef)

  Opt_sd = survival_fn$opt_sd
  opt_ASP = survival_fn$opt_asp
  ASP = sum((opt_ASP-seq(1,0.1,-0.1))^2)
  Surv_time = survival_fn$Surv_time
  return(list(Surv_time = Surv_time, ASP = round(ASP,4), Opt_sd=Opt_sd, opt_ASP = opt_ASP))
}

# Ftn_fit_mdl function
Ftn_fit_mdl = function(reg_model, data, y_var='eventtime', status_var='status'){
  mdl_dat =  data%>% select(-status_var)
  if(toupper(reg_model) == "LIN"){
    mdl = lm(paste0(y_var, " ~."), data = mdl_dat)
    y_hat = calib_survtime(pred_survtime=predict(mdl), min_survtime=min(mdl_dat[,y_var])); BCRF_coef = NULL
  }else if(toupper(reg_model) == "RF"){
    mdl = ranger(paste0(y_var, " ~."), data =  mdl_dat, num.trees=100, seed=1886)
    y_hat = calib_survtime(pred_survtime=mdl$predictions, min_survtime=min(mdl_dat[,y_var])); BCRF_coef = NULL
  }
  return (list(mdl=mdl, y_hat=y_hat))
}    

# Ftn_cal_Surv_time function
Ftn_cal_Surv_time = function(newdata, y_var, status_var, discr_names, BCRF_coef, num_sample, conti_range){
  min_time = min(newdata[,y_var])
  max_time = max(newdata[,y_var])+sd(newdata[,y_var])
  
  if(class(model) == "ranger") {Surv_ftn = RF_Surv_ftn} else
    if(class(model) == "lm") {Surv_ftn = Lin_Surv_ftn} else
      if(class(model) == "xgb.Booster") {Surv_ftn = XGB_Surv_ftn}
  
  conti_idx = which(!(colnames(newdata) %in% c(discr_names, y_var, status_var)))
  discr_idx = which(colnames(newdata) %in% discr_names)
  
  SPR_Surv = Surv_ftn(conti_range=conti_range, fitted_model=model, newdata=newdata,y_var, conti_idx=conti_idx, discr_idx=discr_idx, num_sample = num_sample, sample_seed = 1886, max_time=max_time, min_time=min_time, BCRF_coef=BCRF_coef)
  Surv_time = as.data.frame(t(sapply(1:nrow(newdata), function(x) Check_My_Surv(SPR_Surv[[x]], max_time, seq(1,0.1,-0.1), EoR_time))))
  colnames(Surv_time) = paste0('P', seq(1, 0.1, -0.1))
  rownames(Surv_time) = paste0("Obs",1:nrow(newdata))
  
  return(list(Surv_time = Surv_time, SPR_Surv = SPR_Surv))
}

# Ftn_cal_ASP function
Ftn_cal_ASP = function(SPR_Surv, newdata, y_var, status_var, EoR_time){
  max_time = max(newdata[,y_var])+sd(newdata[,y_var])
  ASP_ = SPR_ASP(SPR_Surv, newdata, y_var, status_var, max_time, EoR_time)
  ASP = sum((ASP_ - seq(1, 0.1, -0.1))^2)
  return(round(ASP,4))
}

