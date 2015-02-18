# Create the Correct Tables for Prothrombin Data
source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

Liver_Data = read.table('./LiverData.csv', sep = ',', header = TRUE)

Table_2 = Liver_Data[,c(1,2,3,4)]

Table_2$treatment[Table_2$time == 0] = -1

Table_2$treatment = as.factor(Table_2$treatment)

levels(Table_2$treatment) = c('null', 'control', 'treatment')

names(Table_2)[3] = 'obs_times'
names(Table_2)[2] = 'obs'

id = as.numeric(levels(as.factor(Liver_Data$id)))

Table_1 = data.frame(id)

id_check <- function(x) {which(Liver_Data$id == x)[1]}

Table_1$cens = 1-Liver_Data$cens[as.numeric(lapply(Table_1$id, id_check))]
Table_1$survival = Liver_Data$survival[as.numeric(lapply(Table_1$id, id_check))]

## Initialization 
delta = 1/365
censored <- function(x) {Table_1$cens[Table_1$id == x]}
survival <- function(x) {Table_1$survival[Table_1$id == x]}

Table_2$cens = as.numeric(lapply(Table_2$id, censored))
Table_2$survival = as.numeric(lapply(Table_2$id, survival))
Table_2$revival =  Table_2$survival - Table_2$obs_times
Table_2$logrev = log(Table_2$revival+delta)
cutoff = (Table_2$revival >= 0.5)
Table_2$spline = (Table_2$revival-0.5)*cutoff

Table_2_uncens = Table_2[Table_2$cens==0,c(1:4,6:9)]
Table_1_uncens = Table_1[Table_1$cens==0,]

### Compute the Covariance Models

Patient <- outer(Table_2_uncens$id, Table_2_uncens$id, "==")  # Patient Indicator Matrix

cov_lambda <- 1.67
  
Patient.ds <- exp(-abs(outer(Table_2_uncens$revival,Table_2_uncens$revival, "-"))/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix

system.time(baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$treatment+Table_2_uncens$survival+Table_2_uncens$revival+Table_2_uncens$logrev, ~Patient + Patient.ds, kernel = 0))

cov_params = baseline_model$sigma
mean_params = baseline_model$beta

system.time(int_baseline_model <- regress(Table_2_uncens$obs~Table_2_uncens$treatment+Table_2_uncens$survival+Table_2_uncens$revival+Table_2_uncens$treatment*Table_2_uncens$logrev, ~Patient + Patient.ds, kernel = 0))

cov_params_int <- int_baseline_model$sigma
beta_int <- int_baseline_model$beta

## For a patient's records, generate Covariate matrix 
## and Sigma given parameters beta and sigma^2

Covariate <- function(data, T) {
	revival = T-data[,1]
	delta = 1/365
	log_rev = log(revival+delta)
	survival = rep(T,length(revival))
  if (length(revival) == 1) {
    return(matrix(c(rep(1,length(revival)),data[,2:3], survival, revival,log_rev), nrow = 1))
  } else {
    return(cbind(rep(1,length(revival)),data[,2:3], survival, revival,log_rev))
  }
}

Covariate_int <- function(data, T) {
  revival = T-data[,1]
  delta = 1/365
  log_rev = log(revival+delta)
  survival = rep(T,length(revival))
  if (length(revival) == 1) {
    return(matrix(c(rep(1,length(revival)),data[,2:3], survival, revival,log_rev,log_rev*data[,2:3]), nrow = 1))
  } else {
    return(cbind(rep(1,length(revival)),data[,2:3], survival, revival,log_rev,log_rev*data[,2:3]))
  }
}

Sigma <- function(data, T, cov_params) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  sigmasq_2 = cov_params[3]
  lambda = 1.67
  obs_times = data[,1]
  k = length(obs_times)
  
  return( sigmasq_0 * diag(k) + sigmasq_1 * outer(rep(1,k),rep(1,k))
          + sigmasq_2 * exp(-abs(outer(obs_times, obs_times,"-"))/lambda))
}