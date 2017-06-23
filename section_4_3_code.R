# Create the Correct Tables for Prothrombin Data

liver <- read.csv("http://www.stat.uchicago.edu/~pmcc/reports/LiverData.csv")

baseline <- liver$time==0
liver$treatment[baseline] <- -1
censored <- !liver$cens		## apparently liver$cens is coded in reverse
revival <- (liver$survival - liver$time)[!censored]
id <- liver$id[!censored];  prothrombin <- liver$prothrombin[!censored];
survival <- liver$survival[!censored];  time <- liver$time[!censored]; 
treat <- liver$treatment[!censored]
treat <- as.factor(treat);  levels(treat) <- c("null", "control", "prednizone")

## Initialization 
fs <- trunc(survival);  fs <- pmin(fs,8)
ft <- trunc(time);  ft <- pmin(ft, 8)
fst <- fs + 9*ft
ymean <- round(tapply(prothrombin, fst, mean), 1)
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
liver$revival = liver$survival - liver$time
r1 <- tapply(revival, id, min)
min2 <- function(x){diff(sort(x))[1]};  r2 <- tapply(revival, id, min2)
min3 <- function(x){diff(sort(x))[2]};  r3 <- tapply(revival, id, min3)
min4 <- function(x){diff(sort(x))[3]};  r4 <- tapply(revival, id, min4)
summary(r1*365);  summary(r2*365); summary(r3*365);  summary(r4*365)