# Functions
log.lik <- function(Table_1, Table_2) {
  internal.llikind <- function(beta) {
    ll.sum = 0
    surv.ids = unique(Table_1$id[Table_1$cens == 0])
    for (id in surv.ids) {
      times = Table_2$obs_times[Table_2$id == id]
      times = times[times != 0]
      surv.time = Table_1$survival[Table_1$id == id]
      treat = Table_1$treat[Table_1$id == id]
      ll.sum = ll.sum + log.lik.ind(beta)(times, treat, surv.time)
    }
    return(-ll.sum)
  }
  return(internal.llikind)
}

log.lik.ind <- function(beta) {
  internal.llikind <- function(times, treat, surv.time) {
    # Compute individual log-likelihood using
    # approximation by days-level constant rate
    
    seq.time = seq(0,surv.time,delta)
    r.time = (seq.time < 1)*log(3) + (seq.time >= 1)*log(1)
    reverse.time = surv.time-seq.time
    invlin.time = reverse.time/(reverse.time+tau)
    X.time = cbind(1, invlin.time, r.time, treat)
    lambda.time = exp(X.time%*%beta)
    
    r.appt.times = (times < 1)*log(3) + (times >= 1)*log(1)
    reverse.appt.times = surv.time - times
    invlin.appt.times = reverse.appt.times/(reverse.appt.times+tau)
    if(length(times) == 0) {
      loglambda.appt.times = 0.0
    } else {
      X.appt.times = cbind(1, invlin.appt.times, r.appt.times, treat)
      loglambda.appt.times = X.appt.times%*%beta
    }
    
    
    return(-sum(lambda.time)*delta + sum(loglambda.appt.times))
  }
  return(internal.llikind)
}

# Create the Correct Tables for Prothrombin Data

Liver_Data = read.table('/Users/walterdempsey/joint_models/revival_models/revival_paper/prednisone code/imputation/LiverData.csv', sep = ',', header = TRUE)

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
Table_1$treat = Liver_Data$treatment[as.numeric(lapply(Table_1$id, id_check))]

##  
delta = 1/365
tau = 0.01
beta = c(0,0,1,0)

temp.llik = log.lik(Table_1, Table_2)

output = optim(beta, temp.llik, lower = c(-10,-10,-10), 
      upper = c(10,10,10) )

output$par
output$value

library(numDeriv)

hess.llik = hessian(temp.llik, output$par)

std.err = sqrt(diag(solve(hess.llik)))

round(cbind(output$par,std.err, output$par/std.err),2)

### Let's plot this

# 2 Years leading up to failure
max.T = 2
seq.time = seq(0,max.T,delta)
revival.time = max.T - seq.time
revival.invlin.time = revival.time/(revival.time+tau)
X.time = cbind(1,revival.invlin.time, 0, 0)

Sigma.par = solve(hess.llik)

log.lambda = output$par[1] + output$par[2] * revival.invlin.time
upperCI.log.lambda = lowerCI.log.lambda = log.lambda*0

for(i in 1:length(log.lambda)) {
  upperCI.log.lambda[i] = log.lambda[i] + 1.96*sqrt(t(c(1,revival.invlin.time[i], 0, 0))%*%Sigma.par%*%c(1,revival.invlin.time[i], 0, 0))
  lowerCI.log.lambda[i] = log.lambda[i] - 1.96*sqrt(t(c(1,revival.invlin.time[i], 0, 0))%*%Sigma.par%*%c(1,revival.invlin.time[i], 0, 0))
}

log.lambda.treat = output$par[1] + output$par[2] * revival.invlin.time + output$par[4]
upperCI.log.lambda.treat = lowerCI.log.lambda.treat = log.lambda.treat*0

for(i in 1:length(log.lambda)) {
  upperCI.log.lambda.treat[i] = log.lambda.treat[i] + 1.96*sqrt(t(c(1,revival.invlin.time[i], 0, 1))%*%Sigma.par%*%c(1,revival.invlin.time[i], 0, 1))
  lowerCI.log.lambda.treat[i] = log.lambda.treat[i] - 1.96*sqrt(t(c(1,revival.invlin.time[i], 0, 1))%*%Sigma.par%*%c(1,revival.invlin.time[i], 0, 1))
}

par(oma = c(2.5,2,0,1) + 0.1,
    mar = c(1,1,0.5,0) + 0.1)
plot(-revival.time, log.lambda, type = "l", ylim = c(0.0,1.5), axes = FALSE)
lines(-revival.time, lowerCI.log.lambda, lty = 2)
lines(-revival.time, upperCI.log.lambda, lty = 2)
axis(side = 1, cex.axis = 0.75)
axis(side = 2, cex.axis = 0.75)
mtext(text = "Time until failure",side = 1, line = 2, cex = 0.75)
mtext(text = "Hazard (Log scale)",side = 2, line = 2, cex = 0.75)

lines(-revival.time, log.lambda.treat, lty = 1, col = "red")

legend(-1.95, 1.5, c("Control", "Treatment"), 
       col = c("black", "red"), 
       lty = c(1, 1), 
       cex = 0.8, bty = "n")


### Empirical bayes estimates

sum(exp(log.lambda[revival.time < 2.0 & revival.time > 1.5]))*delta
sum(exp(log.lambda[revival.time < 1.5 & revival.time > 1.0]))*delta
sum(exp(log.lambda[revival.time < 1.0 & revival.time > 0.5]))*delta
sum(exp(log.lambda[revival.time < 0.5 & revival.time > 0.0]))*delta

