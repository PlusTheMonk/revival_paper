## Survival Model ##
#source('uncens_regress.R')
cov_params = c(1634.790, 0.00, 6244.585)
mean_params = c(67.684083, 2.043645,30.873804,10.617377, 264.775823, -132.034929)
source('../density_components.R')

## Find Survival Distributions ##
divat.data.surv = divat.data[divat.data$Deces == 0 & !is.na(divat.data$creat),]

ids = unique(divat.data.surv$id)
covs = matrix(nrow = 1, ncol =7)
time = cens = rep(0,0)

for (i in ids) {
  record = divat.data.surv[divat.data.surv$id == i, ]
  covs = rbind(covs,
               c(record$centre[1],record$yearGreffe[1],record$rangGreffe[1],
                 record$AgeR[1],record$sexeR[1], record$tailleR[1], 
                 record$poids[1]))
  time = c(time, record$TpsEvtAns_depM12[1])
  cens = c(cens, record$Retour[1] == 0)
}
covs = covs[-1,]

covs[,2] = covs[,2] - min(covs[,2])

age = covs[,4]

library(survival)

exp.fit <- survreg(Surv(time, cens) ~ as.factor(covs[,1]) + 
                    covs[,2],
                   dist="exponential")
summary(exp.fit)

parameters = exp.fit$coefficients; surv_model = 'exponential'

weib.fit <- survreg(Surv(time, cens) ~ as.factor(covs[,1]) + 
                      covs[,2],,
                    dist="weibull")
summary(weib.fit)

# parameters = list('theta' = 1/exp(exp_model$coefficients)); surv_model = 'exponential'
# parameters = list('lambda' = exp(weibull_model$coefficients), 'k' = weibull_model$scale); surv_model = 'weibull'

# Check Exponential Model via Sensitivity using CI
#parameters = list('theta' = 1/exp(log.theta.high)); surv_model = 'exponential'
#parameters = list('theta' = 1/exp(log.theta.low)); surv_model = 'exponential'

## Revival Model ##
formula = ~-1+tps_postM12+sexeR

Data = model.matrix(formula, data = divat.data.surv)

### Imputation Output ####
set.seed(123)

beta_list = beta_cov_list = sigma_list = sigma_cov_list = list()
M = 5; num_imputes = 1

mean_params = baseline_model$beta
cov_params = baseline_model$sigma
surv_model = "exponential"

nondeces.id = unique(divat.data.surv$id)

for (k in 1:M) {
  X_imp = rep(0,0); Y_imp = rep(0,0); id_imp = rep(0,0)
  for (j in nondeces.id[1:75]) {
    if (any(divat.data.surv$Retour[divat.data.surv$id == j] == 0)) {
      pat_table = matrix(Data[divat.data.surv$id == j,],ncol = 2)
      surv.data = c(1,divat.data.surv$AgeR[divat.data.surv$id ==j][1])
      y = divat.data.surv$creat[divat.data.surv$id == j]
      c = divat.data.surv$TpsEvtAns_depM12[divat.data.surv$id == j][1]
      imp_table = impute_Cov(mean_params, cov_params, parameters, surv.data, surv_model, pat_table, y, c, Covariate)
#       temp = score.info(y,imp_table,mean_params, cov_params)
      X_imp = rbind(X_imp,imp_table)
      Y_imp = c(Y_imp,y)
      id_imp = c(id_imp,rep(paste(j,"_",1, sep = ""),length(y)))
    }
  }

  ### Compute the Covariance Models
  
  Patient <- outer(id_imp, id_imp, "==")  # Patient Indicator Matrix
  
  cov_lambda <- 7
  
  Patient.ds <- exp(-abs(outer(X_imp[,4],X_imp[,4], "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix
    
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