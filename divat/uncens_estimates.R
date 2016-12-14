library("beepr")
# Create the Correct Tables for Prothrombin Data
source('http://www.stat.uchicago.edu/~pmcc/courses/regress.R')

divat.data = read.table('./b.long.csv', sep = ';', header = TRUE)


## Initialization
delta = 1/365

divat.data$revival = divat.data$TpsEvtAns_depM12-divat.data$tps_postM12
divat.data$logrev = log(divat.data$revival+delta)
divat.data$invlin = divat.data$revival/(1+divat.data$revival)

divat.data.retour = divat.data[(divat.data$Retour==1)
                               & (!is.na(divat.data$creat))
                              ,]

nrow(divat.data.retour) # Number of measurements
length(unique(divat.data.retour$id)) # Number of patients with graft loss

### Compute the Covariance Models

Patient <- outer(divat.data.retour$id, divat.data.retour$id, "==")  # Patient Indicator Matrix

cov_lambda <- 7.0
nds = -abs(outer(divat.data.retour$revival,divat.data.retour$revival, "-"))
Patient.ds <- exp(nds/cov_lambda) *Patient # Patient Specific Exponential Covariance Matrix
#Patient.ds <- -abs(outer(divat.data.retour$revival,divat.data.retour$revival, "-")) *Patient # Patient Specific Covariance Matrix
ndt = -abs(outer(divat.data.retour$tps_postM12,divat.data.retour$tps_postM12, "-"))

system.time(baseline_model <-
                regress(divat.data.retour$creat~
                            divat.data.retour$TpsEvtAns_depM12+
                            #as.factor(divat.data.retour$rangGreffe)+
                            #divat.data.retour$yearGreffe+
                            #divat.data.retour$AgeR+
                            divat.data.retour$sexeR+
                            #divat.data.retour$tailleR+
                            divat.data.retour$revival+
                            divat.data.retour$invlin+
                            divat.data.retour$logrev,
                        ~Patient + Patient.ds, kernel = 0,
                        pos = c(1,1)))

baseline_model$llik

beep()

summary(baseline_model)

cov_params = baseline_model$sigma
mean_params = baseline_model$beta

## Checking model adequacy
log.s = log(divat.data.retour$revival+delta)
nonlin <- exp(-abs(outer(log.s, log.s, "-"))/cov_lambda) # Non-linear Covariance Matrix

system.time(nonlin_check <-
              regress(divat.data.retour$creat~
                        divat.data.retour$TpsEvtAns_depM12+
                        #as.factor(divat.data.retour$rangGreffe)+
                        #divat.data.retour$yearGreffe+
                        #divat.data.retour$AgeR+
                        divat.data.retour$sexeR+
                        #divat.data.retour$tailleR+
                        divat.data.retour$revival+
                        divat.data.retour$invlin+
                        divat.data.retour$logrev,
                      ~Patient + Patient.ds + nonlin, kernel = 0, 
                      pos = c(0,1,0,1),
                      taper = c(0.3,0.3,0.5)))
beep()
summary(nonlin_check)
nonlin_check$llik - baseline_model$llik

system.time(ndt_check <-
              regress(divat.data.retour$creat~
                        divat.data.retour$TpsEvtAns_depM12+
                        #as.factor(divat.data.retour$rangGreffe)+
                        #divat.data.retour$yearGreffe+
                        #divat.data.retour$AgeR+
                        divat.data.retour$sexeR+
                        #divat.data.retour$tailleR+
                        divat.data.retour$revival+
                        divat.data.retour$invlin+
                      divat.data.retour$logrev,
                      ~Patient + Patient.ds + ndt, kernel = 0, start = c(cov_params,0)))
beep()
summary(ndt_check)

2*(ndt_check$llik - baseline_model$llik)

## For a patient's records, generate Covariate matrix
## and Sigma given parameters beta and sigma^2

Covariate <- function(data, T) {
	revival = T-data[,2]
	delta = 1/365
	log_rev = log(revival+delta)
  invlin = revival/(revival+1)
	survival = rep(T,length(revival))
  if (length(revival) == 1) {
    return(matrix(c(rep(1,length(revival)),data[,2], survival, revival,invlin, log_rev), nrow = 1))
  } else {
    return(cbind(rep(1,length(revival)),survival, data[,2],revival,invlin, log_rev))
  }
}

Sigma <- function(data, T, cov_params) {
  sigmasq_0 = cov_params[1]
  sigmasq_1 = cov_params[2]
  sigmasq_2 = cov_params[3]
  lambda = 7
  obs_times = data[,2]
  k = length(obs_times)

  return( sigmasq_0 * diag(k) + sigmasq_1 * outer(rep(1,k),rep(1,k))
          + sigmasq_2 * exp(-abs(outer(obs_times, obs_times,"-"))/lambda))
}

