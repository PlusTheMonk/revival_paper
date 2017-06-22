### Functions and Source Code
## Set up
liver.data = read.table("../imputation/LiverData.csv", sep = ",", header = TRUE)

# load packages 'JM' and 'lattice'
library("JM")
library("lattice")

liver.data$sq.time = liver.data$time^2

### Joint Model

##linear mixed effects model
fitLME <- lme(prothrombin ~ time , random = ~time | id, data = liver.data)

##Cox model surv(time, death)~1
patients = unique(liver.data$id)
liver.data.id = liver.data[1:length(patients),]
for (i in 1:length(patients)) {
  patient.data = (liver.data[liver.data$id==patients[i],])
  liver.data.id[i,] = patient.data[1,]
}

# Patient 463 with issue
liver.data.id$survival[liver.data.id$id == 463] = max(liver.data$time[liver.data$id == 463]) +0.00001

fitSURV <- coxph(Surv(survival, cens) ~ 1, data = liver.data.id, x = TRUE)

##method "weibull-AFT-GH"
fit.JM <- jointModel(fitLME, fitSURV, timeVar = "time")
summary(fit.JM)

## Get Patient Records
patient.id = 402;
delta = 0.001
times = seq(max(liver.data$time[liver.data$id == patient.id])+delta,7.5,delta)
plot(liver.data$time[liver.data$id == patient.id], 
     liver.data$prothrombin[liver.data$id == patient.id])

ND.59 <- liver.data[liver.data$id %in% patient.id, ]
ND.59$survival= min(times)
predSurv.59 <- survfitJM(fit.JM, newdata = ND.59, idVar = "id", 
                          last.time = "survival",
                          survTimes = times)

ND.69 <- ND.59; ND.69[8,2] = 69
predSurv.69 <- survfitJM(fit.JM, newdata = ND.69, idVar = "id", 
                          last.time = "survival",
                          survTimes = times)

ND.79 <- ND.59; ND.79[8,2] = 79
predSurv.79 <- survfitJM(fit.JM, newdata = ND.79, idVar = "id", 
                          last.time = "survival",
                          survTimes = times)

surv.59 = c(1,predSurv.59$summaries$`402`[,2])
surv.69 = c(1,predSurv.69$summaries$`402`[,2])
surv.79 = c(1,predSurv.79$summaries$`402`[,2])

hazard.59 = (1-surv.59[2:length(surv.59)]/surv.59[1:(length(surv.59)-1)])/delta
hazard.69 = (1-surv.69[2:length(surv.69)]/surv.69[1:(length(surv.69)-1)])/delta
hazard.79 = (1-surv.79[2:length(surv.79)]/surv.79[1:(length(surv.79)-1)])/delta

lambda = 0.2
png("srem_liver_prediction.png", width = 8,height = 4, units = "in", res = 300)
par(mar = c(3,3,1,1)+0.1)
plot(times[-1],hazard.59, type = "l", axes = FALSE, ylab = "", xlab = "", ylim = c(0.05,0.14),xlim = c(5,7))
lines(times[-1],hazard.69)
lines(times[-1],hazard.79)
axis(side=1, at = seq(5,7,0.5)); axis(side=2)

segments(min(times[-1]), lambda, max(times[-1]), lty = 2)

mtext("Conditional hazard", side = 2,line = 2)
mtext("Time", side = 1,line = 2)
text(6.75,0.135,"59", cex = 0.9)
text(6.75,0.123,"69", cex = 0.9)
text(6.75,0.105,"79", cex = 0.9)

dev.off()


### Prediction function
predict.JM <- function(patient.data,s,t) {
  prediction.data = patient.data[patient.data$tps_postM12< t,]
  prediction.data$TpsEvtAns_depM12 = s
  prediction.data$id= 1
  predSurv <- survfitJM(fit.JM, newdata = prediction.data, idVar = "id", 
                            last.time = "TpsEvtAns_depM12",
                            survTimes = t+s)
  return(predSurv$summaries$`1`[2])
}
