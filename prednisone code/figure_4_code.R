###### Example ############
#load("/Users/walterdempsey/Documents/stat/research/joint_models/ideas/dual-alignment-project/code/prednisone/Liver.RData")

liver <- read.csv("./imputation/LiverData.csv")
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

### Revival Uncensored

for(j in 1:length(rev)) {
	abs_diff = abs(revival-rev[j])
	dmin_time = tapply(abs_diff,id,min)
	response = vector(length = length(dmin_time))

	for (i in 1:num_id) {
		obs = which(abs_diff[id == id_levels[i]]  == dmin_time[i])
		response[i] = mean((prothrombin[id == id_levels[i]])[obs])
	}
	for (k in 1:2) {
		val_rev[j,k] = sum( dnorm(dmin_time/h) * response * (id_treat == k) * (rev[j] <= max_time)) / sum(dnorm(dmin_time/h)* (rev[j] <= max_time) * (id_treat == k)*(rev[j] <= max_time))
	}
}

#### Revival Uncensored ####

cens_revival <- (liver$survival-liver$time)[censored]
cens_forward <- liver$time[censored]

cens_id <- liver$id[censored];  cens_prothrombin <- liver$prothrombin[censored];
cens_survival <- liver$survival[censored];  cens_time <- liver$time[censored];
cens_treat <- liver$treatment[censored]
cens_treat <- as.factor(cens_treat);  levels(cens_treat) <- c("control", "prednizone")
cens_max_time = tapply(cens_revival, cens_id, max)

cens_val_rev = matrix(nrow = length(rev), ncol = 2)
cens_val_forw = matrix(nrow = length(forw), ncol = 2)
cens_id_levels = levels(as.factor(cens_id))
cens_num_id = length(levels(as.factor(cens_id)))
cens_id_treat = tapply(as.numeric(cens_treat),cens_id,min)

for(j in 1:length(rev)) {
  abs_diff = abs(cens_revival-rev[j])
  dmin_time = tapply(abs_diff,cens_id,min)
  response = vector(length = length(dmin_time))
  
  for (i in 1:cens_num_id) {
    obs = which(abs_diff[cens_id == cens_id_levels[i]]  == dmin_time[i])
    response[i] = mean((cens_prothrombin[cens_id == cens_id_levels[i]])[obs])
  }
  for (k in 1:2) {
    cens_val_rev[j,k] = sum( dnorm(dmin_time/h) * response * (cens_id_treat == k) * (rev[j] <= cens_max_time)) / sum(dnorm(dmin_time/h)* (rev[j] <= cens_max_time) * (cens_id_treat == k) * (rev[j] <= cens_max_time))
  }
}

#### Revival Plot ####
png("/Users/walterdempsey/Documents/stat/research/joint_models/revival_models/rev_mean_traj.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(oma = c(2,0,0,2) + 0.1,
          mar = c(1,1,0.5,1) + 0.1)

plot(-rev,cens_val_rev[,1], col = "black", type = "l", lty = 2, ylim = c(50, 100), axes = FALSE, ylab = "Prothrombin Index", xlab = "Time Since Recruitment")
axis(side = 1, cex.axis = 0.75)
axis(side = 4, cex.axis = 0.75)
lines(-rev,cens_val_rev[,2], col = "red", type = "l", lty = 2)

# legend(-2,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), bty = "n", cex = 0.7)
legend(-4,60, c("Complete","Censored"), col = c("black","black"), lty = c(1,2), bty = "n", cex = 0.7)
mtext("Time Until Failure", side = 1,line = 2, cex = 0.75)

lines(-rev,val_rev[,1], col = "black", type = "l", lty = 1)
lines(-rev,val_rev[,2], col = "red", type = "l", lty = 1)

# legend(-4,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), cex = 0.7, bty = "n")

dev.off()

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
		val_forw[j,k] = sum( dnorm(dmin_time/h) * response * (id_treat == k) * (forw[j] <= max_time)) / sum(dnorm(dmin_time/h)* (id_treat == k)*(forw[j] <= max_time))
	}
}

#### Forward Censored ####

for(j in 1:length(forw)) {
  abs_diff = abs(cens_forward-forw[j])
  dmin_time = tapply(abs_diff,cens_id,min)
  response = vector(length = length(dmin_time))
  
  for (i in 1:cens_num_id) {
    obs = which(abs_diff[cens_id == cens_id_levels[i]]  == dmin_time[i])
    response[i] = mean((cens_prothrombin[cens_id == cens_id_levels[i]])[obs])
  }
  for (k in 1:2) {
    cens_val_forw[j,k] = sum( dnorm(dmin_time/h) * response * (cens_id_treat == k) * (forw[j] <= cens_max_time)) / sum(dnorm(dmin_time/h) * (cens_id_treat == k) * (forw[j] <= cens_max_time))
  }
}

#### Forward Plot ####

png("/Users/walterdempsey/Documents/stat/research/joint_models/revival_models/forw_mean_traj.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(oma = c(2,1,0,1) + 0.1,
          mar = c(1,1,0.5,0) + 0.1)

plot(forw,cens_val_forw[,1], col = "black", type = "l", lty = 2, ylim = c(50, 100), axes = FALSE, ylab = "Prothrombin Index", xlab = "Time Since Recruitment")
axis(side = 1, cex.axis = 0.75)
axis(side = 2, cex.axis = 0.75)
lines(forw,cens_val_forw[,2], col = "red", type = "l", lty = 2)

legend(1,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), bty = "n", cex = 0.7)
mtext("Time Since Recruitment", side = 1,line = 2, cex = 0.75)

lines(forw,val_forw[,1], col = "black")
lines(forw,val_forw[,2], col = "red")

# legend(-4,60, c("Control","Prednisone"), col = c("black","red"), lty = c(1,1), cex = 0.7, bty = "n")

dev.off()