### Functions and Source Code
source('uncens_regress.R')
source('density_components.R')

Covariate <- function(data, T) {
  revival = T-data[,1]
  delta = 1/365
  log_rev = log(revival+delta)
  invlin = revival/(revival+1)
  survival = rep(T,length(revival))
  if (length(revival) == 1) {
    return(matrix(c(rep(1,length(revival)),survival, data[,2], revival,invlin, log_rev), nrow = 1))
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

pred <- function( mean_params, cov_params, pat_table,y,Cov) {
  # Provides a function of the survival time, t,
  # for the likelihood Y | T
  
  pred2 <- function(t) {
    k = dim(pat_table)[1]
    X = Cov(pat_table,t)
    S = Sigma(pat_table, t, cov_params)
    mu = X%*%mean_params
    return( -t(y - mu)%*% solve(S) %*% (y - mu)/2 )
  }
  
  return(pred2)	
}

## Set up
formula = ~-1+tps_postM12+sexeR

Data = model.matrix(formula, data = divat.data)

## Get Patient Records
id.patient = 48 # 72
plot(divat.data$tps_postM12[divat.data$id == j], divat.data$creat[divat.data$id == j])
mean_params = baseline_model$beta
cov_params = baseline_model$sigma

png("/Users/walterdempsey/Documents/stat/research/joint_models/revival_models/divat/divat_prediction.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(mfrow = c(1,2),
          oma = c(2,1,0,0) + 0.1,
          mar = c(0,1,1,1) + 0.1)

pat_table = matrix(Data[divat.data$id == id.patient,], ncol = 2)
y = divat.data$creat[divat.data$id == id.patient]
t = divat.data$tps_postM12[divat.data$id == id.patient]
c = divat.data$TpsEvtAns_depM12[divat.data$id == id.patient][1]
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
eval_T = seq(c+0.001,15,0.01)
probs = cond_dens(eval_T)
max_T = round(c+9,0)
min_T = c+0.1
p.minp = probs-probs[length(probs)]
plot(eval_T[eval_T<max_T & eval_T > min_T],(p.minp)[eval_T<max_T  & eval_T > min_T], type = "l", ylim = c(-20,5), xlim = c(round(min_T,0),max_T), axes = FALSE, ylab = "", xlab = "")
axis(side = 1,at = seq(6,max_T,1), cex.axis = 0.6); axis(side = 2, labels = FALSE)
text(7,-0.14,y[length(y)], cex = 0.5)

y2 = y
y2[length(y2)] = 310
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y2, Covariate))
probs = cond_dens(eval_T)
p.minp = probs-probs[length(probs)]
lines(eval_T[eval_T<max_T & eval_T > min_T],(p.minp)[eval_T<max_T & eval_T > min_T]) 
text(6.5,3.85,y2[length(y2)], cex = 0.5)

y3 = y
y3[length(y3)] = 110
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y3, Covariate))
probs = cond_dens(eval_T)
p.minp = probs-probs[length(probs)]
lines(eval_T[eval_T<max_T & eval_T > min_T],(p.minp)[eval_T<max_T & eval_T > min_T]) 
text(8.25,-4,y3[length(y3)], cex = 0.5)

### Hazard Functions : Using the exponential survival function
lambda = 0.02623713
pat_table = matrix(Data[divat.data$id == id.patient,], ncol = 2)
y = divat.data$creat[divat.data$id == id.patient]
max_T = 15
eval_T = seq(min_T,100,0.005)
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
probs = cond_dens(eval_T)
total_probs_unnorm = exp(probs-lambda*eval_T)
total_probs = total_probs_unnorm/sum(total_probs_unnorm)
hazard = total_probs/(1-cumsum(total_probs))

plot(eval_T[eval_T<max_T],(hazard)[eval_T<max_T], type = "l", ylim = c(0,0.0025), xlim = c(round(c,0),max_T), axes = FALSE)
axis(side = 1,at = seq(6,max_T,1), cex.axis = 0.6); axis(side = 2, labels = FALSE)
text(14,0.0005,y[length(y)], cex = 0.5)

segments(min(eval_T),lambda*0.005,x1 = max_T, lty = 2)

y[length(y)] = 110
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
probs = cond_dens(eval_T)
total_probs_unnorm = exp(probs-lambda*eval_T)
total_probs = total_probs_unnorm/sum(total_probs_unnorm)
hazard = total_probs/(1-cumsum(total_probs))
lines(eval_T[eval_T<max_T],hazard[eval_T<max_T])
text(14,0.00025,y[length(y)], cex = 0.5)

y[length(y)] = 310
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
probs = cond_dens(eval_T)
total_probs_unnorm = exp(probs-lambda*eval_T)
total_probs = total_probs_unnorm/sum(total_probs_unnorm)
hazard = total_probs/(1-cumsum(total_probs))
lines(eval_T[eval_T<max_T],hazard[eval_T<max_T])
text(14,0.0009,y[length(y)], cex = 0.5)


dev.off()