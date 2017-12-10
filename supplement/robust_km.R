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

### Standard Weibull Fit 
library('survival')
expfit = survreg(Surv(survival, !cens) ~ 1, Table_1, dist='weibull',
                 scale=1)

weibfit = survreg(Surv(survival, !cens) ~ 1, Table_1, dist='weibull')

exp.lambda = exp(expfit$coefficients)

weib.lambda = exp(weibfit$coefficients)
weib.gamma = 1/weibfit$scale

### Tail only fit

times = unique(Table_1$survival[Table_1$cens == 0])
maxT = max(Table_1$survival)

hat.kappa = 1/(log(maxT) - mean(log(times)))
hat.lambda = maxT/length(times)^(1/hat.kappa)

hat.exp.lambda = maxT/length(times)
### 

seqt = seq(0,15,0.1)

weib.surv = exp(-(seqt/weib.lambda)^weib.gamma)
exp.surv = exp(-(seqt/exp.lambda))

kap_meier_adj <- function(lambda,kappa) {
  temp_f <- function(t) {
    obs = times < t
    maxT = max(Table_1$survival)
    logtot = 0
    for(s in times[obs]) {
      d = sum(Table_1$survival[Table_1$cens == 0] == s)
      r = sum(Table_1$survival[Table_1$cens == 0] > s) +
          sum(Table_1$survival[Table_1$cens == 1] >= s)
      logtot = logtot + log(r) - log(r+d)
    }
    if (t > maxT) {
      addtot = exp((maxT/lambda)^kappa - (t/lambda)^kappa)
    } else {addtot = 1}
    return(exp(logtot)*addtot)
  }
  return(temp_f)
}

weib.km.surv = Vectorize(kap_meier_adj(hat.lambda, hat.kappa))(seqt)
exp.km.surv = Vectorize(kap_meier_adj(hat.exp.lambda, 1))(seqt)

par(mar = c(3.5,2,1,2)+0.1)
plot(seqt, weib.surv, type= "l", axes = FALSE, ylim = c(0,1), xlab = "", ylab = "")
lines(seqt,exp.surv, col = "red")
lines(seqt,weib.km.surv, lty = 2)
#lines(seqt[seqt > maxT],exp.km.surv[seqt>maxT], lty = 2, col = "red")
axis(side = 1, cex.axis = 0.75); axis(side = 2, cex.axis = 0.75)
mtext(text="Time (in years)", side=1, line=2, cex = 0.75)
legend(10,0.8,c("Exponential", "Weibull", expression(rho==0)), col = c("red","black","black"), lty= c(1,1,2),cex = 0.75, bty = "n")

### Second example

epsilon = 0.001
group_1 = c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
group_0 = c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35)

cens_1 = rep(0,21)
cens_0 = c(0,0,0,1,0,1,0,1,1,0,0,1,1,1,0,0,1,1,1,1,1)

times = unique(group_0[cens_0 == 0])
maxT = max(group_0)

hat.kappa = 1/(log(maxT) - mean(log(times)))
hat.lambda = maxT/length(times)^(1/hat.kappa)

hat.exp.lambda = maxT/length(times)

## Fits 
expfit = survreg(Surv(group_0, !cens_0) ~ 1, dist='weibull',
                 scale=1)

weibfit = survreg(Surv(group_0, !cens_0) ~ 1, dist='weibull')

exp.lambda = exp(expfit$coefficients)

weib.lambda = exp(weibfit$coefficients)
weib.gamma = 1/weibfit$scale

## Plot for Gehan

seqt = seq(0,50,0.1)

weib.surv = exp(-(seqt/weib.lambda)^weib.gamma)
exp.surv = exp(-(seqt/exp.lambda))

kap_meier_adj_gehan <- function(lambda,kappa) {
  temp_f <- function(t) {
    obs = times < t
    maxT = max(group_0)
    logtot = 0
    for(s in times[obs]) {
      d = sum(group_0[cens_0 == 0] == s)
      r = sum(group_0[cens_0 == 0] > s) +
        sum(group_0[cens_0 == 1] >= s)
      logtot = logtot + log(r) - log(r+d)
    }
    if (t > maxT) {
      addtot = exp((maxT/lambda)^kappa - (t/lambda)^kappa)
    } else {addtot = 1}
    return(exp(logtot)*addtot)
  }
  return(temp_f)
}

weib.km.surv.gehan = Vectorize(kap_meier_adj_gehan(hat.lambda, hat.kappa))(seqt)
exp.km.surv.gehan = Vectorize(kap_meier_adj_gehan(hat.exp.lambda, 1))(seqt)

par(mar = c(3.5,2,1,2)+0.1)
plot(seqt, weib.surv, type= "l", axes = FALSE, ylim = c(0,1), xlab = "")
lines(seqt,exp.surv, col = "red")
lines(seqt,weib.km.surv.gehan, lty = 2)
#lines(seqt[seqt > maxT],exp.km.surv.gehan[seqt > maxT], lty = 2, col = "red")
axis(side = 1, cex.axis = 0.75); axis(side = 2, cex.axis = 0.75)
mtext(text="Time (in weeks)", side=1, line=2, cex = 0.75)
legend(30,0.9,c("Exponential", "Weibull", expression(rho==0)), col = c("red","black","black"), lty= c(1,1,2),cex = 0.75, bty = "n")
