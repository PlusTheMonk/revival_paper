observations = rep(0,0)
n = 200
set.seed(1)

for (i in 1:n) {
  T = rexp(1, rate = 1/5)
  
  t = seq(0,floor(T))
  id = rep(i, length(t))
  s = T-t
  kept_apt = runif(length(t)) <= 5/(5+t)
  
  Sigma = (diag(length(t))+1+exp(- abs(outer(t,t,"-"))/5))/2
  mu = 10 + 10*s/(10+s)
  
  Y = mu+t(chol(Sigma))%*%rnorm(length(t))
  
  observations = rbind(observations,cbind(id,Y,t,s,kept_apt))
  
}

Y = observations[,2]
t = observations[,3]
s = observations[,4]
kept_apt = observations[,5]

# png("/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/simulated_health.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(mfrow = c(1,2),
          oma = c(2,1,0,2) + 0.1,
          mar = c(0,1,1,1) + 0.1)

plot(t[as.logical(kept_apt)],Y[as.logical(kept_apt)], axes = FALSE, cex = 0.75, xlim = c(0,20), ylim = c(6,20))
lines(lowess(t[as.logical(kept_apt)],Y[as.logical(kept_apt)]), col = "red", lwd = 4)
axis(side = 1, cex.axis = 0.75); axis(side = 2, cex.axis = 0.75)


plot(-s[as.logical(kept_apt)],Y[as.logical(kept_apt)], axes = FALSE, cex = 0.75, xlim = c(-20,0), ylim = c(6,20))
lines(lowess(-s[as.logical(kept_apt)],Y[as.logical(kept_apt)]), col = "red", lwd = 4)
axis(side = 1, cex.axis = 0.75); axis(side = 4, cex.axis = 0.75)

# dev.off()

### Conditional Standard Deviation Calculations

cond_dens <- function(data) {
  cond_dens_t <- function(T) {
    kept_apt = data[,5]
    obs_times = data[kept_apt==1,3]
    obs = data[kept_apt==1,2]
    Sigma = (diag(length(obs_times))+1+exp(- abs(outer(obs_times,obs_times,"-"))/5))/2    
    rev = T - obs_times
    mu = 10 + 10*rev/(10+rev)
    return ( exp(-t(obs-mu)%*%solve(Sigma,obs-mu)/2)*dexp(T,rate = 1/5) )
  }
  return(cond_dens_t)
}


std_dev = pred_err = exp_pred_err = rep(0,0)

for (i in 1:n) {
  data = matrix(observations[observations[,1] == i,], ncol = 5)
  if (sum(data[,5]) >= 2) {
    kept_apt = data[,5]
    eval_T = seq(max(data[kept_apt==1,3]),max(data[kept_apt==1,3])+20,0.01)
    prob = unlist(lapply(eval_T,cond_dens(data)))
    norm_prob = prob/sum(prob)
    mu = sum(norm_prob*eval_T)
    mu_exp = max(data[kept_apt==1,3])+5
    std_dev = c(std_dev,sqrt(sum((eval_T-mu)^2*norm_prob)))
    survival = data[1,3]+data[1,4]
    pred_err = c(pred_err,(survival-mu)^2)
    exp_pred_err = c(exp_pred_err,(survival-mu_exp)^2)
  }
}


RMSE = sqrt(mean(pred_err))
exp_RMSE = sqrt(mean(exp_pred_err))
Avg_StdDev = mean(std_dev)

print(RMSE)
print(exp_RMSE)
print(Avg_StdDev)