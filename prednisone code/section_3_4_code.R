# Create the Correct Tables for Prothrombin Data

Liver_Data = read.table('./imputation/LiverData.csv', sep = ',', header = TRUE)

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

## Initialization 
delta = 1/365
censored <- function(x) {Table_1$cens[Table_1$id == x]}
survival <- function(x) {Table_1$survival[Table_1$id == x]}

Table_2$cens = as.numeric(lapply(Table_2$id, censored))
Table_2$survival = as.numeric(lapply(Table_2$id, survival))
Table_2$revival =  Table_2$survival - Table_2$obs_times

# computation of Table 1 (mean prothrombin values)
fs <- trunc(Table_2$survival[Table_2$cens == 0]);  fs <- pmin(fs,8)
ft <- trunc(Table_2$obs_times[Table_2$cens == 0]);  ft <- pmin(ft, 8)
fst <- fs + 9*ft
ymean <- round(tapply(Table_2$obs[Table_2$cens == 0], fst, mean), 1)
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
Table_2_uncens = Table_2[Table_2$cens == 0,]
id_uncens = unique(Table_2_uncens$id)

revinterval <- function(i) {
  rev = Table_2_uncens$revival[Table_2_uncens$id == i]     
  rev_ints = rev[length(rev)]
  if (length(rev) > 1) {
    rev_ints = c(-rev[2:length(rev)] + rev[1:(length(rev)-1)],rev_ints)  
  }
  return(rev_ints)
}

rev_ints = lapply(id_uncens,revinterval)

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

