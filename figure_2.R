# Create the Correct Tables for Prothrombin Data

divat.data = read.table('./b.long.csv', sep = ';', header = TRUE)

deces.ind = unique(divat.data$id[divat.data$Deces == 1])

divat.data$revival = divat.data$TpsEvtAns_depM12-divat.data$tps_postM12

rev.kernelmean <- function(obs, h = 0.3) {
    revfn <- function(s) {

    weight = dnorm((divat.data$revival[obs]-s)/h)

    return(sum(weight*divat.data$creat[obs])/sum(weight))
}
    return(revfn)
}

forw.kernelmean <- function(obs, h = 0.3) {
    forwfn <-  function(t) {    
    
    sub.obs = obs & (divat.data$TpsEvtAns_depM12 > t)

    weight = dnorm((divat.data$tps_postM12[sub.obs]-t)/h)

    return(sum(weight*divat.data$creat[sub.obs])/sum(weight))

}
    return(forwfn)
}

seq.time = seq(0.1,4,0.1)

deces.obs = (divat.data$Deces == 1) & (!is.na(divat.data$creat))
retour.obs = (divat.data$Retour == 1) & (!is.na(divat.data$creat))
cens.obs = (divat.data$Retour == 0) & (divat.data$Deces == 0) & (!is.na(divat.data$creat))


rev.creat.deces = Vectorize(rev.kernelmean(deces.obs))(seq.time)
forw.creat.deces = Vectorize(forw.kernelmean(deces.obs))(seq.time)

rev.creat.retour = Vectorize(rev.kernelmean(retour.obs))(seq.time)
forw.creat.retour = Vectorize(forw.kernelmean(retour.obs))(seq.time)

rev.creat.cens = Vectorize(rev.kernelmean(cens.obs))(seq.time)
forw.creat.cens = Vectorize(forw.kernelmean(cens.obs))(seq.time)


# Forward plot
png("forw_divat_traj.png")
par(mar = c(3,3,1,1)+0.1, cex = 0.75)
plot(seq.time,forw.creat.retour, type = "l", ylim = c(120,400), axes = FALSE)
lines(seq.time, forw.creat.cens, col = "red")
lines(seq.time, forw.creat.deces, lty = 2)
axis(side = 2); axis(side = 1)
mtext("Time Since 1-year post transplantation", side = 1, cex = 0.75, line=2)
mtext(expression(paste("Serum creatinine (",mu,"mol/L)")), side = 2, cex = 0.75, line=2)
dev.off()
# Revival plot
png("rev_divat_traj.png")
par(mar = c(3,2,1,3)+0.1, cex = 0.75)
plot(-seq.time,rev.creat.retour, type = "l", ylim = c(120,400), axes = FALSE)
lines(-seq.time, rev.creat.cens, col = "red")
lines(-seq.time, rev.creat.deces, lty = 2)
axis(side = 4); axis(side = 1)
mtext("Time Until Failure", side = 1, cex = 0.75, line=2)
legend(-4,400,c("Graft Loss", "Deceased", "Censored"), col = c("black", "black", "red"), lty = c(1,2,1),bty ="n", cex = 1.5)
dev.off()