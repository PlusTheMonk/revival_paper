# Create the Correct Tables for Prothrombin Data

divat.data = read.table('./b.long.csv', sep = ';', header = TRUE)

## Initialization 
delta = 1/365
deces.ind = unique(divat.data$id[divat.data$Retour == 1])

surv <- function(id) {
    divat.data$TpsEvtAns_depM12[divat.data$id == id][1]
}

deces.surv = unlist(lapply(deces.ind,surv))

# computation of Table 1 (mean prothrombin values)
obs = ((divat.data$Retour == 1) & (!is.na(divat.data$creat)))
fs <- trunc(divat.data$TpsEvtAns_depM12[obs]);  fs <- pmin(fs,8)
ft <- trunc(divat.data$tps_postM12[obs]);  ft <- pmin(ft, 8)
fst <- fs + 9*ft
ymean <- round(tapply(divat.data$creat[obs], fst, mean), 1)
M <- N <- matrix(0, 9, 9); row <- as.numeric(gl(9, 1, 81));  col <- as.numeric(gl(9, 9, 81))
colnames(M) <- colnames(N) <- c("t0-1", "t1-2", "t2-3", "t3-4", "t4-5", "t5-6", "t6-7", "t7-8", "t8+")
rownames(M) <- rownames(N) <- c("T0-1", "T1-2", "T2-3", "T3-4", "T4-5", "T5-6", "T6-7", "T7-8", "T8+")
M[row >= col] <- ymean;  N[row >= col] <- table(fst)

print(M)  # Table 1
print(N)

#####
ymean <- as.vector(M)
rtime <- row - col
w <- as.numeric(row >= col) 
fit0 <- lm(ymean~1, weights=w);  rss0 <- sum(w*fit0$resid^2)
fitforw <- lm(ymean~as.factor(col), weights=w);  rssforw <- sum(w*fitforw$resid^2)
fitrev <- lm(ymean~as.factor(rtime), weights=w);  rssrev <- sum(w*fitrev$resid^2)

rss0-rssforw
rss0-rssrev

### ANOVA
fitnull <- lm(ymean~as.factor(col) + as.factor(rtime) + as.factor(row), weights=w);  rssnull <- sum(w*fitnull$resid^2)
fitcr <- lm(ymean~as.factor(col) + as.factor(row), weights=w);  rsscr <- sum(w*fitcr$resid^2); dfcr = 7
fitcd <- lm(ymean~as.factor(col) + as.factor(rtime), weights=w);  rsscd <- sum(w*fitcd$resid^2); dfcd = 7
fitrd <- lm(ymean~as.factor(rtime) + as.factor(row), weights=w);  rssrd <- sum(w*fitrd$resid^2); dfrd = 7
fit_res <- lm(ymean~as.factor(col)*as.factor(row), weights=w);  rss_res <- sum(w*fit_res$resid^2); df_res = 21

round(rbind(c(rsscr - rssnull, dfcr, (rsscr-rssnull)/dfcr),
c(rssrd-rssnull, dfrd, (rssrd-rssnull)/dfrd),
c(rsscd-rssnull, dfcd, (rsscd-rssnull)/dfcd),
c(rssnull-rss_res, df_res, rssnull/df_res)),1)

### Interval Check : Section 2.2
divat.data$revival = divat.data$TpsEvtAns_depM12-divat.data$tps_postM12
divat.data.uncens = divat.data[divat.data$Retour == 1,]

revinterval <- function(i) {
  rev = divat.data.uncens$revival[divat.data.uncens$id == i]     
  rev_ints = rev[length(rev)]
  if (length(rev) > 1) {
    rev_ints = c(-rev[2:length(rev)] + rev[1:(length(rev)-1)],rev_ints)  
  }
  return(rev_ints)
}

rev_ints = lapply(deces.ind,revinterval)

final = second_final = third_final = rep(0,0)

for(i in 1:length(rev_ints)){
  x = rev_ints[[i]]
  final = c(final,x[length(x)])
  
  if(length(x) > 1 ) {
    second_final = c(second_final,x[length(x)-1])   
  }
  if(length(x) > 2 ) {
    third_final = c(third_final,x[length(x) - 2])
  }
}


c(median(final), median(second_final), median(third_final))*365
c(mean(final), mean(second_final), mean(third_final))*365

