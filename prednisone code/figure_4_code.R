###### Example ############

setwd("/Users/walterdempsey/Documents/stat/research/joint_models/dual-alignment-project/code/prednisone")
load("Liver.RData")

source("http://www.stat.uchicago.edu/~pmcc/courses/regress.R")

liver <- read.csv("LiverData.csv")
censored <- !liver$cens

### Non Paramateric UnCensored Curves ###

#### Revival Plot ####
revival <- (liver$survival-liver$time)[!censored]
forward <- liver$time[!censored]

id <- liver$id[!censored];  prothrombin <- liver$prothrombin[!censored];
survival <- liver$survival[!censored];  time <- liver$time[!censored];
treat <- liver$treatment[!censored]
treat <- as.factor(treat);  levels(treat) <- c("control", "prednizone")
max_time = tapply(revival, id, max)


h = 0.5
forw = rev = seq(0,5,0.2)
val_rev = matrix(nrow = length(rev), ncol = 2)
val_forw = matrix(nrow = length(forw), ncol = 2)
id_levels = levels(as.factor(id))
num_id = length(levels(as.factor(id)))
id_treat = tapply(as.numeric(treat),id,min)

for(j in 1:length(rev)) {
	abs_diff = abs(revival-rev[j])
	dmin_time = tapply(abs_diff,id,min)
	response = vector(length = length(dmin_time))

	for (i in 1:num_id) {
		obs = which(abs_diff[id == id_levels[i]]  == dmin_time[i])
		response[i] = mean((prothrombin[id == id_levels[i]])[obs])
	}
	for (k in 1:2) {
		val_rev[j,k] = sum( dnorm(dmin_time/h) * response * (id_treat == k) * (t[j] <= max_time)) / sum(dnorm(dmin_time/h)* (t[j] <= max_time) * (id_treat == k))
	}
}


#### Forward Uncensored ####
max_time = tapply(forward, id, max)

for(j in 1:length(forw)) {
	abs_diff = abs(forward-forw[j])
	dmin_time = tapply(abs_diff,id,min)
	response = vector(length = length(dmin_time))

	for (i in 1:num_id) {
		obs = which(abs_diff[id == id_levels[i]]  == dmin_time[i])
		response[i] = mean((prothrombin[id == id_levels[i]])[obs])
	}
	for (k in 1:2) {
		val_forw[j,k] = sum( dnorm(dmin_time/h) * response * (id_treat == k) * (t[j] <= max_time)) / sum(dnorm(dmin_time/h)* (t[j] <= max_time) * (id_treat == k))
	}
}

png("/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/prot_mean_trajectories.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(mfrow = c(1,2),
          oma = c(2,1,0,2) + 0.1,
          mar = c(1,1,0.5,1) + 0.1)

plot(forw,val_forw[,1], col = "black", type = "l", lty = 1, ylim = c(50, 90), axes = FALSE, ylab = "Prothrombin Index", xlab = "Time Since Recruitment")
axis(side = 1, cex.axis = 0.75)
axis(side = 2, cex.axis = 0.75)
lines(forw,val_forw[,2], col = "red", type = "l", lty = 2)

legend(1,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), bty = "n", cex = 0.7)
mtext("Time Since Recruitment", side = 1,line = 2, cex = 0.75)


plot(-rev,val_rev[,1], col = "black", type = "l", lty = 1, ylim = c(50, 90), axes = FALSE, ylab = "Prothrombin Index", xlab = "Time Until Failure")
axis(side = 1, cex.axis = 0.75)
axis(side = 4, cex.axis = 0.75)
lines(-rev,val_rev[,2], col = "red", type = "l", lty = 2)
mtext("Time Until Failure", side = 1,line = 2, cex = 0.75)

# legend(-4,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), cex = 0.7, bty = "n")

dev.off()