# Functions
log.lik.selfcorrecting <- function(Table_1, Table_2) {
  internal.llikind <- function(beta) {
    ll.sum = 0
    surv.ids = unique(Table_1$id[Table_1$cens == 0])
    for (id in surv.ids) {
      times = Table_2$obs_times[Table_2$id == id]
      times = times[times != 0]
      surv.time = Table_1$survival[Table_1$id == id]
      treat = Table_1$treat[Table_1$id == id]
      ll.sum = ll.sum + log.lik.ind.selfcorrecting(beta)(times, treat, surv.time)
    }
    return(-ll.sum)
  }
  return(internal.llikind)
}

log.lik.ind.selfcorrecting <- function(beta) {
  internal.llikind <- function(times, treat, surv.time) {
    # Compute individual log-likelihood using
    # approximation by days-level constant rate
    
    seq.time = seq(0,surv.time,delta)
    # r.time = (seq.time < 1)*log(3) + (seq.time >= 1)*log(1)
    withinyear.count <- function (t) {
      sum(
        (times[times != 0] >= floor(t)) & 
          (times[times != 0] < t )
      )
    }
    
    N.time = unlist(lapply(seq.time, withinyear.count))
    self.correct = (seq.time - floor(seq.time) - 3*N.time) * (seq.time <= 1) + 
      (seq.time - floor(seq.time) - 1*N.time) * (seq.time > 1)
    r.time = (seq.time < 1)*log(3) + (seq.time >= 1)*log(1)
    reverse.time = surv.time-seq.time
    invlin.time = reverse.time/(reverse.time+tau)
    X.time = cbind(1, invlin.time, r.time, self.correct, treat)
    lambda.time = exp(X.time%*%beta)
    
    N.appt.times = unlist(lapply(times, withinyear.count))
    self.correct.appt.times = (times - floor(times) - 3*N.appt.times) * (times <= 1) + 
      (times - floor(times) - 1*N.appt.times) * (times > 1)
    r.appt.times = (times < 1)*log(3) + (times >= 1)*log(1)
    reverse.appt.times = surv.time - times
    invlin.appt.times = reverse.appt.times/(reverse.appt.times+tau)
    if(length(times) == 0) {
      loglambda.appt.times = 0.0
    } else {
      X.appt.times = cbind(1, invlin.appt.times, r.appt.times, self.correct.appt.times, treat)
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
beta = c(4,-4,0.6, -0.1,0)

temp.llik = log.lik.selfcorrecting(Table_1, Table_2)

temp.llik(beta)

output.sc = optim(beta, temp.llik, lower = rep(-5,4),
               upper = rep(5,4), control=list(trace=TRUE) )

output.sc$par
output.sc$value

# 4.09588457 -3.94953416  1.03718761  0.24036306 -0.03160661
# 336.1587

library(numDeriv)

hess.llik = hessian(temp.llik, output.sc$par)

std.err.sc = sqrt(diag(solve(hess.llik)))

round(cbind(output.sc$par,std.err.sc, output.sc$par/std.err.sc),2)

# [1,]  4.10       0.14  29.74
# [2,] -3.95       0.14 -28.62
# [3,]  1.04       0.05  18.96
# [4,]  0.24       0.02  12.83
# [5,] -0.03       0.05  -0.58

### Let's plot this

# 2 Years leading up to failure
max.T = 2
seq.time = seq(0,max.T,delta)
revival.time = max.T - seq.time
revival.invlin.time = revival.time/(revival.time+tau)
times = Table_2$obs_times[Table_2$id == 402]
times = times[times!= 0]
treat = Table_1$treat[Table_1$id == 402]

beta.sc = output.sc$par
beta = output$par
eval_T = seq(max(times)+0.001,15,0.01)
llik.apptschedule = llik.apptschedule.sc = eval_T*0
for(i in 1:length(eval_T)) {
  llik.apptschedule.sc[i] = log.lik.ind.selfcorrecting(beta)(times, treat, eval_T[i])  
  llik.apptschedule[i] = log.lik.ind(beta)(times, treat, eval_T[i])  
}

plot(eval_T[eval_T<7],(llik.apptschedule.sc-llik.apptschedule.sc[length(llik.apptschedule.sc)])[eval_T<7], type = "l", ylim = c(6,18), xlim = c(5,7), axes = FALSE, ylab = "", xlab = "")
lines(eval_T[eval_T<7],(llik.apptschedule-llik.apptschedule[length(llik.apptschedule)])[eval_T<7])