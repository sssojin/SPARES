# SPARES

You need to install packages using `install.packages`:
- survival
- randomForestSRC
- ranger
- xgboost
- dplyr
- e1071
- parallel
- doParallel
- GeneralizedHyperbolic
- dplyr

The SPARES model consists of two modules: SPARES_Imputation and SPARES_ASP.
Please refer to the file `SPARES_run.Rmd`.

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

# SPARES_ASP

# Input
# - dat : output from the SPARES_Imputation
# - EoR_time : time of end of research

# Output
# - surv_time : survival probability matrix for each observation
# - asp : asp for the model

# -----------------------------------------------------------------







