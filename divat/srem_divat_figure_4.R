### Functions and Source Code
## Set up
divat.data = read.table("../b.long.csv", sep = ";", header = TRUE)

divat.data = divat.data[-which(is.na(divat.data$creat)),]

# load packages 'JM' and 'lattice'
library("JM")
library("lattice")

### Joint Model

##linear mixed effects model
fitLME <- lme(creat ~ tps_postM12 , random = ~ tps_postM12 | id, data = divat.data)

##Cox model surv(time, death)~1
patients = unique(divat.data$id)
divat.data.id = divat.data[1:length(patients),]
for (i in 1:length(patients)) {
  patient.data = (divat.data[divat.data$id==patients[i],])
  divat.data.id[i,] = patient.data[1,]
}

fitSURV <- coxph(Surv(TpsEvtAns_depM12, Retour) ~ 1, data = divat.data.id, x = TRUE)

##method "weibull-AFT-GH"
fit.JM <- jointModel(fitLME, fitSURV, timeVar = "tps_postM12")
summary(fit.JM)

## Get Patient Records
patient.id = 48;
times = seq(6.2260,15.0,0.001)
plot(divat.data$tps_postM12[divat.data$id == patient.id], 
     divat.data$creat[divat.data$id == patient.id])

ND.210 <- divat.data[divat.data$id %in% patient.id, ]
predSurv.210 <- survfitJM(fit.JM, newdata = ND.210, idVar = "id", 
                          last.time = "TpsEvtAns_depM12",
                          survTimes = times)

ND.110 <- ND.210; ND.110[7,3] = 110
predSurv.110 <- survfitJM(fit.JM, newdata = ND.110, idVar = "id", 
                          last.time = "TpsEvtAns_depM12",
                          survTimes = times)

ND.310 <- ND.210; ND.310[7,3] = 310
predSurv.310 <- survfitJM(fit.JM, newdata = ND.310, idVar = "id", 
                          last.time = "TpsEvtAns_depM12",
                          survTimes = times)

surv.110 = c(1,predSurv.110$summaries$`48`[,2])
surv.210 = c(1,predSurv.210$summaries$`48`[,2])
surv.310 = c(1,predSurv.310$summaries$`48`[,2])

hazard.210 = 1-surv.210[2:length(surv.210)]/surv.210[1:(length(surv.210)-1)]
hazard.310 = 1-surv.310[2:length(surv.310)]/surv.310[1:(length(surv.310)-1)]
hazard.110 = 1-surv.110[2:length(surv.110)]/surv.110[1:(length(surv.110)-1)]


png("/Users/walterdempsey/Documents/stat/research/joint_models/revival_models/divat/srem_divat_prediction.png", width = 6.5,height = 3, units = "in", res = 300)
par(mar = c(1,1,1,1)+0.1)
plot(times[-1],hazard.210, type = "l", axes = FALSE, ylab = "", xlab = "", ylim = c(0,0.0004),xlim = c(6,15))
lines(times[-1],hazard.110)
lines(times[-1],hazard.310)
axis(side=1, at = seq(6,15,1),cex.axis = 0.75); axis(side=2,labels =  FALSE, cex.axis = 0.75)

segments(min(times[-1]), lambda*0.001, max(times[-1]), lty = 2)

mtext("Conditional hazard", side = 2,line = 1, cex = 0.75)
mtext("Time", side = 1,line = 2, cex = 0.75)
text(8,0.000085,"210", cex = 0.7)
text(8,0.00016,"310", cex = 0.7)
text(8,0.000045,"110", cex = 0.7)

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
