library(survival); library(randomForestSRC);
library(dplyr); library(gridExtra); library(parallel); library(doParallel)

# 3 Scenarios Running (1 PH & 2 nonPHs)

# Scenario 1) Cox Proportional Hazards model --------------------------------------------------------

# ----------------------------------------------------
# Main function for generating simulation data
# ----------------------------------------------------
sim_run_PH = function(distn, para, covs, cov.betas){
  if(distn == "exponential"){
    survT = apply(covs, 1, function(x) -log(runif(1)) / (para[1] * exp(sum(cov.betas * x))))
  }else if(distn == "weibull"){
    survT = apply(covs, 1, function(x) (-log(runif(1)) / (para[1] * exp(sum(cov.betas * x))))^(1/para[2]))
  }else if(distn == "gompertz"){
    survT = apply(covs, 1, function(x) {1/para[2] * log(1 - para[2]*log(runif(1))/(para[1] * exp(sum(cov.betas * x))))})
  }
  return(data.frame(eventtime = survT, covs))
}


## Survival error sum calculation function - coxPH
cal_PH_ES = function(sim.dat, PH.predict, distn = NULL, para = NULL, betas, predictors){
  Time = PH.predict$time
  rec_len = diff(Time)
  if(is.null(distn)){
    St = sim.dat$baseline$survivor[sim.dat$baseline$time %in% Time]
    sim.dat = sim.dat$data
  }
  
  cal_ES = function(ind){
    elp = exp(sum(sim.dat[ind, predictors] * betas))
    
    if(is.null(distn)){
      true_St = St[-length(Time)]^elp
    }else if(distn == "exponential"){
      true_St = exp(-para[1] * Time * elp)[-length(Time)]
    }else if(distn == "weibull"){
      true_St = exp(-para[1]*Time^para[2] * elp)[-length(Time)]
    }else if(distn == "gompertz"){
      true_St = (1-pgompertz(Time[-length(Time)], shape = para[2], rate = para[1]))^elp
    }
    fitted_St = PH.predict$surv[,ind][-length(Time)]
    
    abs_err = sum(abs((true_St - fitted_St) * rec_len))
    square_err = sum((true_St - fitted_St)^2 * rec_len)
    return(c(abs_err, square_err))
  }
  ES = t(sapply(1:nrow(sim.dat), cal_ES))
  colnames(ES) = c("Absolute_Err", "Square_Err")
  return(ES)
}

## Survival error sum calculation function - RSF
cal_RSF_ES = function(sim.dat, RSF.predict, distn = NULL, para = NULL, betas, predictors){
  Time = RSF.predict$time.interest
  rec_len = diff(Time)
  if(is.null(distn)){
    St = sim.dat$baseline$survivor[sim.dat$baseline$time %in% Time]
    sim.dat = sim.dat$data
  }
  
  cal_ES = function(ind){
    elp = exp(sum(sim.dat[ind, predictors] * betas))
    
    if(is.null(distn)){
      true_St = St[-length(Time)]^elp
    }else if(distn == "exponential"){
      true_St = exp(-para[1] * Time[-length(Time)] * elp)
    }else if(distn == "weibull"){
      true_St = exp(-para[1]*Time^para[2] * elp)[-length(Time)]
    }else if(distn == "gompertz"){
      true_St = (1-pgompertz(Time[-length(Time)], shape = para[2], rate = para[1]))^elp
    }
    fitted_St = RSF.predict$survival[ind,-length(Time)]
    
    abs_err = sum(abs((true_St - fitted_St) * rec_len))
    square_err = sum((true_St - fitted_St)^2 * rec_len)
    return(c(abs_err, square_err))
  }
  ES = t(sapply(1:nrow(sim.dat), cal_ES))
  colnames(ES) = c('Absolute_Err', "Square_Err")
  return(ES)
}


# Scenario 2) Cox Proportional Hazards model with time-varying covariates ---------------------------

## Survival times generation function
expo_T = function(x, xc.ind, elp.ind, scale, betas, beta.t){
  1/(beta.t*x[xc.ind]) * log(1 + beta.t*x[xc.ind]*(-log(runif(1))) / (scale*exp(sum(x[elp.ind] * betas))))
}

gomp_T = function(x, xc.ind, elp.ind, scale, shape, betas, beta.t){
  1/(beta.t*x[xc.ind] + shape) * log(1 + (beta.t*x[xc.ind] + shape) * (-log(runif(1))) / (scale*exp(sum(x[elp.ind] * betas))))
}

## Cumulative hazard function to calculate survival probability
expo_HtX = function(t, x, xc.ind, elp.ind, scale, betas, beta.t){
  scale*exp(sum(x[elp.ind] * betas)) / (beta.t*x[xc.ind]) * (exp(beta.t*x[xc.ind]*t)-1)
}

gomp_HtX = function(t, x, xc.ind, elp.ind, scale, shape, betas, beta.t){
  scale*exp(sum(x[elp.ind] * betas)) / (beta.t*x[xc.ind] + shape) * (exp((beta.t*x[xc.ind] + shape)*t)-1)
}


# ----------------------------------------------------
# Main function for generating simulation data
# ----------------------------------------------------
sim_run_nonPH1 = function(distn, para, covs, xc.ind, cov.betas, beta.t, ...){
  elp.ind = dplyr::setdiff(1:ncol(covs), xc.ind)
  if(distn == "exponential"){
    survT = apply(covs, 1, function(x) expo_T(x, xc.ind, elp.ind, para[1], cov.betas, beta.t))
  }else if(distn == "gompertz"){
    survT = apply(covs, 1, function(x) gomp_T(x, xc.ind, elp.ind, para[1], para[2], cov.betas, beta.t))
  }
  return(data.frame(eventtime = survT, covs))
}


## Survival error sum calculation function - coxPH
cal_PH_ES2 = function(sim.dat, PH.predict, distn, para, betas, beta.t, k = NULL, predictors){
  Time = PH.predict$time
  rec_len = diff(Time)
  p = length(predictors)
  
  cal_ES = function(ind){
    cal.mu = sim.dat[ind, predictors]
    
    if(distn == "exponential"){
      Ht = expo_HtX(Time, unlist(cal.mu), para[1], betas, beta.t, k, p)
      true_St = exp(-Ht[-length(Time)])
    }else if(distn == "gompertz"){
      Ht = gomp_HtX(Time, unlist(cal.mu), para[1], para[2], betas, beta.t, k, p)
      true_St = exp(-Ht[-length(Time)])
    }
    fitted_St = PH.predict$surv[,ind][-length(Time)]
    
    abs_err = sum(abs((true_St - fitted_St) * rec_len))
    square_err = sum((true_St - fitted_St)^2 * rec_len)
    return(c(abs_err, square_err))
  }
  ES = t(sapply(1:nrow(sim.dat), cal_ES))
  colnames(ES) = c("Absolute_Err", "Square_Err")
  return(ES)
}

## Survival error sum calculation function - RSF
cal_RSF_ES2 = function(sim.dat, RSF.predict, distn, para, betas, beta.t, k = NULL, predictors){
  Time = RSF.predict$time.interest
  rec_len = diff(Time)
  p = length(predictors)
  
  cal_ES = function(ind){
    cal.mu = sim.dat[ind, predictors]
    
    if(distn == "exponential"){
      Ht = expo_HtX(Time, unlist(cal.mu), para[1], betas, beta.t, k, p)
      true_St = exp(-Ht[-length(Time)])
    }else if(distn == "gompertz"){
      Ht = gomp_HtX(Time, unlist(cal.mu), para[1], para[2], betas, beta.t, k, p)
      true_St = exp(-Ht[-length(Time)])
    }
    fitted_St = RSF.predict$survival[ind,-length(Time)]
    
    abs_err = sum(abs((true_St - fitted_St) * rec_len))
    square_err = sum((true_St - fitted_St)^2 * rec_len)
    return(c(abs_err, square_err))
  }
  ES = t(sapply(1:nrow(sim.dat), cal_ES))
  colnames(ES) = c('Absolute_Err', "Square_Err")
  return(ES)
}


# Scenario 3) Cox Proportional Hazards model with time-varying and interacted covariates ------------

# ----------------------------------------------------
# Main function for generating simulation data
# ----------------------------------------------------
sim_run_nonPH2 = function(distn, para, covs, discr.ind, xc.ind, cov.betas, beta.t, ...){
  if(distn == "exponential"){
    PH_T = function(x) {-log(runif(1)) / (para[1] * exp(sum(cov.betas * x)))}
    nonPH_T = function(x) expo_T(x, xc.ind, elp.ind = 1:ncol(covs), para[1], cov.betas, beta.t)
  }else if(distn == "gompertz"){
    PH_T = function(x) {1/para[2] * log(1 - para[2]*log(runif(1))/(para[1] * exp(sum(cov.betas * x))))}
    nonPH_T = function(x) gomp_T(x, xc.ind, elp.ind = 1:ncol(covs), para[1], para[2], cov.betas, beta.t) 
  }
  
  survT = apply(covs, 1, function(x) if(x[discr.ind] == 0) PH_T(x) else nonPH_T(x))
  
  return(data.frame(eventtime = survT, covs))
}


# Scenario 4) Cox Proportional Hazards model with switching treatment -------------------------------

# ----------------------------------------------------
# Main function for generating simulation data
# ----------------------------------------------------
sim_run_nonPH3 = function(distn, para, covs, cov.betas, beta.t, ts, xc = NULL, ...){
  p = ncol(covs)
  
  if(is.null(xc)){
    if(distn == "exponential"){
      survT = apply(covs, 1, function(x) expo_Ts(x, para[1], cov.betas, beta.t, ts))
    }else if(distn == "weibull"){
      survT = apply(covs, 1, function(x) weib_Ts(x, para[1], para[2], cov.betas, beta.t, ts))
    }else if(distn == "gompertz"){
      survT = apply(covs, 1, function(x) gomp_Ts(x, para[1], para[2], cov.betas, beta.t, ts))
    }
  }else{
    covs = cbind(covs, xc = xc)
    if(distn == "exponential"){
      survT = apply(covs, 1, function(x) expo_Ts(x[1:p], para[1], cov.betas, beta.t*x[p+1], ts))
    }else if(distn == "weibull"){
      survT = apply(covs, 1, function(x) weib_Ts(x[1:p], para[1], para[2], cov.betas, beta.t*x[p+1], ts))
    }else if(distn == "gompertz"){
      survT = apply(covs, 1, function(x) gomp_Ts(x[1:p], para[1], para[2], cov.betas, beta.t*x[p+1], ts))
    }
    covs = covs[,-ncol(covs)]
  }
  
  return(data.frame(eventtime = survT, covs))
}

## Survival times generation function
expo_Ts = function(x, scale, betas, beta.t, ts){
  u = runif(1)
  if(-log(u) < scale*exp(sum(x * betas))*ts){
    -log(u) / (scale*exp(sum(x * betas)))
  }else{
    (-log(u)-scale*exp(sum(x * betas))*ts+scale*exp(sum(x * betas) + beta.t)*ts) / (scale*exp(sum(x * betas) + beta.t))
  }
}

weib_Ts = function(x, scale, shape, betas, beta.t, ts){
  u = runif(1)
  if(-log(u) < scale*exp(sum(x * betas))*ts^shape){
    (-log(u) / (scale*exp(sum(x * betas))))^(1/shape)
  }else{
    (1/exp(beta.t) * (-log(u)/(scale*exp(sum(x * betas))) + (exp(beta.t)-1) * ts^shape))^(1/shape)
  }
}

gomp_Ts = function(x, scale, shape, betas, beta.t, ts){
  u = runif(1)
  if(-log(u) < scale*exp(sum(x * betas))/shape * (exp(shape*ts)-1)){
    1/shape * log(1 + shape*(-log(u)) / (scale*exp(sum(x * betas))))
  }else{
    1/shape * log(shape*(-log(u)) / (scale*exp(sum(x * betas) + beta.t)) - (exp(shape*ts)-1-exp(beta.t + shape*ts)) / exp(beta.t))
  }
}


# Add noise to simulation data ----------------------------------------------------------------------

add_noise = function(sim.org.dat, y.var, sample_seed){
  n = nrow(sim.org.dat)
  set.seed(sample_seed); shuff = sample(1:n)
  set.seed(NULL); noise = rgamma(n, 0.5, 6)
  
  sim.org.dat[,y.var] = sapply(1:n, function(idx){
    if(idx %in% shuff[1:round(n/2)]){
      return(sim.org.dat[idx,y.var]+noise[idx])
    }else{
      Negative = sim.org.dat[idx,y.var]-noise[idx]
      if(Negative >= 0) return(Negative) else return(sim.org.dat[idx,y.var])
    }
  })
  return(sim.org.dat)
}



# Censoring adjustment ------------------------------------------------------------------------------

# Censoring adjusted function
censor_adj = function(N, censorship, sim.dat, y_col, maxt, EoR.prop){
  # the number of samples
  n.EoR = as.integer(EoR.prop * N)
  n.censor = as.integer(censorship * N)
  
  # Separate by maxt
  Lmaxt = which(sim.dat[,y_col] < maxt)
  Umaxt = which(sim.dat[,y_col] >= maxt)
  
  if(length(Umaxt) < n.censor) return(NULL)
  
  # Sample by each case
  d.idx = sample(Lmaxt, N-n.censor)
  RC.idx = sample(dplyr::setdiff(Lmaxt, d.idx), n.censor-n.EoR)
  EoR.idx = sample(Umaxt, n.EoR)
  
  # Adjust eventtime
  censor_y = sapply(sim.dat[RC.idx,y_col], function(t) runif(1, 0, t))
  sim.dat[RC.idx,y_col] = censor_y
  sim.dat[EoR.idx,y_col] = maxt
  
  return(rbind(
    cbind(sim.dat[d.idx,], status = 1, idx = d.idx),
    cbind(sim.dat[c(EoR.idx, RC.idx),], status = 0, idx = c(EoR.idx, RC.idx))
  ))
}

# Censoring adjusted function with weighted probability
censor_adj_wgt = function(censorship, sim.dat, y_col, maxt, EoR.prop){
  n.censor = as.integer(censorship * nrow(sim.dat))
  maxt_idx = which(sim.dat[,y_col] > maxt)
  
  # EoR indicator (0: not censored, 1: random censoring, 2: EoR censoring)
  EoR = rep(0, nrow(sim.dat)); EoR[maxt_idx] = 2
  
  if(censorship == EoR.prop){
    # Censorship = EoR.prop case
    if(abs(mean(EoR == 2)-EoR.prop) > 0.015){
      return(NULL)
    }else{
      sim.dat[maxt_idx,y_col] = maxt
      censor.ind = rep(1, nrow(sim.dat)); censor.ind[maxt_idx] = 0
      EoR = ifelse(censor.ind == 0, 2, 0)
      return(data.frame(sim.dat, status = censor.ind, EoR = EoR))
    }
  }
  
  # Censorship /= EoR.prop case
  if((length(maxt_idx) >= n.censor) | (abs(mean(EoR == 2)-EoR.prop) > 0.02) | length(maxt_idx) == 0){
    return(NULL)
  }else{
    sim.dat[maxt_idx,y_col] = maxt
    
    if(n.censor-length(maxt_idx) > 0){
      weight_seq = seq(0, maxt, length = 11)
      censor_cand = dplyr::setdiff(1:nrow(sim.dat), maxt_idx)
      interval_ind = sapply(1:10, function(x){
        sapply(sim.dat[censor_cand,y_col], function(y){
          return((weight_seq[x] <= y) & y < weight_seq[x+1])
        })
      })
      censor_weight = apply(interval_ind, 1, which)
      
      censor_idx = sample(censor_cand, n.censor-length(maxt_idx), prob = censor_weight)
      censor_y = sapply(sim.dat[censor_idx,y_col], function(t) runif(1, 0, t))
      
      sim.dat[censor_idx,y_col] = censor_y
      EoR[censor_idx] = 1
      maxt_idx = c(maxt_idx, censor_idx)
    }
    
    censor.ind = rep(1, nrow(sim.dat)); censor.ind[maxt_idx] = 0
    return(data.frame(sim.dat, status = censor.ind, EoR = EoR))
  }
}


# For Results ---------------------------------------------------------------------------------------

# Data set for survival function comparison
surv_comp = function(coxPH, RSF, cal.mu){
  PH.predict = summary(survfit(coxPH, newdata = cal.mu))
  PH = data.frame(Time = PH.predict$time, PH.surv = PH.predict$surv)
  
  RSF.predict = predict(RSF, newdata = cal.mu)
  RSF = data.frame(Time = RSF.predict$time.interest, RSF.surv = RSF.predict$survival[1,])
  
  return(list(PH = PH, RSF = RSF))
}


# Percentile Graph
pct_graph = function(comp, censor, p){
  ggplot() +
    geom_line(data = comp[[3]], aes(x = Time, y = True.surv, lty = "True"), size = 1) +
    geom_line(data = comp[[1]], aes(x = Time, y = PH.surv, lty = "coxPH"), size = 1) +
    geom_line(data = comp[[2]], aes(x = Time, y = RSF.surv, lty = "RSF"), size = 1) +
    scale_linetype_manual(values = c("True" = "solid", "coxPH" = "twodash", "RSF" = "dotted")) +
    labs(y = "Survival probability", linetype = "Model",
         title = paste0(censor * 100, "% Censorship"),
         subtitle = paste0(p * 100, "% percentile subject")) +
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.text = element_text(size = 12)
    )
}