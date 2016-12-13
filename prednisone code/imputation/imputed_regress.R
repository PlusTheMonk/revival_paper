## Survival Model ##
source('uncens_regress.R')
source('density_components.R')

library('survival')
exp_model = survreg(Surv(Table_1$survival, !Table_1$cens)~1, dist = 'exponential')
weibull_model = survreg(Surv(Table_1$survival, !Table_1$cens)~1, dist = 'weibull')

parameters = list('theta' = 1/exp(exp_model$coefficients)); surv_model = 'exponential'
# parameters = list('lambda' = exp(weibull_model$coefficients), 'k' = weibull_model$scale); surv_model = 'weibull'

## Revival Model ##
formula = ~-1+obs_times+treatment

Data = model.matrix(formula, data = Table_2)

id = Table_2$id

Data = Data[,c(1,3:4)]

id_rev = Table_2$id
id_surv = Table_1$id

V = Table_1$survival
Y = Table_2$obs
cens = Table_2$cens

### Imputation Output ####
set.seed(123)

beta_list = beta_cov_list = sigma_list = sigma_cov_list = list()
M = 5; num_imputes = 1

mean_params = baseline_model$beta
cov_params = baseline_model$sigma

for (k in 1:M) {
  X_imp = rep(0,0); Y_imp = rep(0,0); id_imp = rep(0,0)
  for (j in 1:488) {
    if (Table_1$cens[Table_1$id == j] == 1) {
      pat_table = matrix(Data[id == j,], ncol = 3)
      y = Y[id == j]
      c = Table_1$survival[Table_1$id == j]
      imp_table = impute_Cov(mean_params, cov_params, parameters, surv_model, pat_table, y, c, Covariate)
      X_imp = rbind(X_imp,imp_table)
      Y_imp = c(Y_imp,y)
      id_imp = c(id_imp,rep(paste(j,"_",1, sep = ""),length(y)))
    }
  }

  ### Compute the Covariance Models
  
  Patient <- outer(id_imp, id_imp, "==")  # Patient Indicator Matrix
  
  cov_lambda <- 1.67
  
  Patient.ds <- exp(-abs(outer(X_imp[,5],X_imp[,5], "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix
    
  imp_model <- regress(Y_imp~X_imp, ~Patient + Patient.ds, kernel = 0)
  
  beta_list[[k]] = imp_model$beta
  beta_cov_list[[k]] = imp_model$beta.cov
  sigma_list[[k]] = imp_model$sigma
  sigma_cov_list[[k]] = imp_model$sigma.cov
}

### Combining Estimates
avg_beta = avg_sigma = 0
avg_beta_cov = avg_sigma_cov = 0
for (k in 1:M) {
  avg_beta = avg_beta + beta_list[[k]]/M
  avg_beta_cov = avg_beta_cov + beta_cov_list[[k]]/M
  avg_sigma = avg_sigma + sigma_list[[k]]/M
  avg_sigma_cov = avg_sigma_cov + sigma_cov_list[[k]]/M
}

betw_beta_cov = betw_sigma_cov = 0

for (k in 1:M) {
  x = (beta_list[[k]] - avg_beta)
  betw_beta_cov = betw_beta_cov + matrix(x%*%t(x), ncol = length(x))/(M-1)
  y = (sigma_list[[k]] - avg_sigma)
  betw_sigma_cov = betw_sigma_cov + matrix(y%*%t(y), ncol = length(y))/(M-1)
}
  
V_beta = avg_beta_cov + (1+1/M)*betw_beta_cov

V_sigma = avg_sigma_cov + (1+1/M)*betw_sigma_cov

beta_imp_se = sqrt(diag(V_beta)); sigma_imp_se = sqrt(diag(V_sigma))

T_stat = (avg_beta-baseline_model$beta)/sqrt(diag(avg_beta_cov))

Imp_beta = cbind(avg_beta,beta_imp_se,avg_beta/beta_imp_se)

Regress_beta = cbind(mean_params,baseline_model$beta.se,mean_params/baseline_model$beta.se)

Imp_sigma = cbind(avg_sigma,sigma_imp_se)

Regress_sigma = cbind(cov_params,sqrt(diag(baseline_model$sigma.cov)))


### Table 4/5 or 8/9 : Depending on Choice of Weibull or Exponential

# (avg_beta-baseline_model$beta)/sqrt(V_beta+baseline_model$beta.se^2)

round(cbind(Imp_beta,Regress_beta,T_stat),2)
round(cbind(Imp_sigma,Regress_sigma),2)

### Interaction Model
set.seed(123)

intbeta_list = intbeta_cov_list = intsigma_list = intsigma_cov_list = list();
M = 5

mean_params = beta_int
cov_params = cov_params_int

for (k in 1:M) {
  intX_imp = rep(0,0); Y_imp = rep(0,0); id_imp = rep(0,0)
  for (j in 1:488) {
    if (Table_1$cens[Table_1$id == j] == 1) {
      pat_table = matrix(Data[id == j,], ncol = 3)
      y = Y[id == j]
      c = Table_1$survival[Table_1$id == j]
      
      imp_table = impute_Cov(mean_params, cov_params, parameters, surv_model, pat_table, y, c, Covariate_int)
      intX_imp = rbind(intX_imp,imp_table)
      Y_imp = c(Y_imp,y)
      id_imp = c(id_imp,rep(paste(j,"_",1, sep = ""),length(y)))  
    }
  }
  
  
  Patient <- outer(id_imp, id_imp, "==")  # Patient Indicator Matrix
  
  cov_lambda <- 1.67
  
  Patient.ds <- exp(-abs(outer(intX_imp[,5],intX_imp[,5], "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix
  
  intimp_model <- regress(Y_imp~intX_imp, ~Patient + Patient.ds, kernel = 0)

  intbeta_list[[k]] = intimp_model$beta
  intbeta_cov_list[[k]] = intimp_model$beta.cov
  intsigma_list[[k]] = intimp_model$sigma
  intsigma_cov_list[[k]] = intimp_model$sigma.cov
}

### Combining Estimates
avg_beta_int = avg_sigma_int = 0
avg_beta_cov_int = avg_sigma_cov_int = 0
for (k in 1:M) {
  avg_beta_int = avg_beta_int + intbeta_list[[k]]/M
  avg_beta_cov_int = avg_beta_cov_int + intbeta_cov_list[[k]]/M
  avg_sigma_int = avg_sigma_int + intsigma_list[[k]]/M
  avg_sigma_cov_int = avg_sigma_cov_int + intsigma_cov_list[[k]]/M
}

betw_beta_cov_int = betw_sigma_cov_int = 0

for (k in 1:M) {
  x = (intbeta_list[[k]] - avg_beta_int)
  betw_beta_cov_int = betw_beta_cov_int + matrix(x%*%t(x), ncol = length(x))/(M-1)
  y = (intsigma_list[[k]] - avg_sigma_int)
  betw_sigma_cov_int = betw_sigma_cov_int + matrix(y%*%t(y), ncol = length(y))/(M-1)
}

V_beta_int = avg_beta_cov_int + (1+1/M)*betw_beta_cov_int

V_sigma_int = avg_sigma_cov_int + (1+1/M)*betw_sigma_cov_int

beta_imp_se_int = sqrt(diag(V_beta_int)); sigma_imp_se_int = sqrt(diag(V_sigma_int))

cons_T_stat_int = (avg_beta_int-int_baseline_model$beta)/sqrt(diag(V_beta_int)+diag(int_baseline_model$beta.cov))
T_stat_int = (avg_beta_int-int_baseline_model$beta)/sqrt(diag(avg_beta_cov_int))

Imp_beta_int = cbind(avg_beta_int,beta_imp_se_int,avg_beta_int/beta_imp_se_int)

Regress_beta_int = cbind(int_baseline_model$beta,int_baseline_model$beta.se,int_baseline_model$beta/int_baseline_model$beta.se)

Imp_sigma_int = cbind(avg_sigma_int,sigma_imp_se_int)

Regress_sigma_int = cbind(int_baseline_model$sigma,sqrt(diag(int_baseline_model$sigma.cov)))

### Table 4/5 or 8/9 : Depending on Choice of Weibull or Exponential

# (avg_beta-baseline_model$beta)/sqrt(V_beta+baseline_model$beta.se^2)

round(cbind(Imp_beta_int,Regress_beta_int,T_stat_int),2)
round(cbind(Imp_sigma_int,Regress_sigma_int),2)
