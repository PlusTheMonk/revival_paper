### Functions and Source Code
source('/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/Code/prednisone code/imputation/uncens_regress.R')
source('/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/Code/prednisone code/imputation/density_components.R')

pred <- function( mean_params, cov_params, pat_table,y,Cov) {
  # Provides a function of the survival time, t,
  # for the likelihood Y | T
  
  pred2 <- function(t) {
    k = dim(pat_table)[1]
    X = Cov(pat_table,t)
    S = Sigma(pat_table, t, cov_params)
    mu = X%*%mean_params
    return( -t(y - mu)%*% solve(S) %*% (y - mu)/2 )
  }
  
  return(pred2)	
}

## Set up
formula = ~-1+obs_times+treatment

Data = model.matrix(formula, data = Table_2)

id = Table_2$id

Data = Data[,c(1,3:4)]

id_rev = Table_2$id
id_surv = Table_1$id

V = Table_1$survival
Y = Table_2$obs
cens = Table_2$cens


## Get Patient Records
j = 402
mean_params = baseline_model$beta
cov_params = baseline_model$sigma

# png("/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/prediction.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(mfrow = c(1,2),
          oma = c(2,1,0,0) + 0.1,
          mar = c(0,1,1,1) + 0.1)

pat_table = matrix(Data[id == j,], ncol = 3)
y = Y[id == j]
c = Table_1$survival[Table_1$id == j]
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)

plot(eval_T[eval_T<7],(probs-probs[length(probs)])[eval_T<7], type = "l", ylim = c(0,2), xlim = c(5,7), axes = FALSE, ylab = "", xlab = "")
axis(side = 1,at = seq(5,7,0.5), cex.axis = 0.75); axis(side = 2, labels = FALSE)
text(eval_T[1],max(probs-probs[length(probs)])*1.05,"59", cex = 0.5)

y[8] = 69
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)
lines(eval_T[eval_T<7],(probs-probs[length(probs)])[eval_T<7]) 
text(eval_T[1],max(probs-probs[length(probs)])*1.05,"69", cex = 0.5)

y[8] = 79
cond_dens = Vectorize(pred(mean_params, cov_params, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)
lines(eval_T[eval_T<7],(probs-probs[length(probs)])[eval_T<7]) 
text(eval_T[1],max(probs-probs[length(probs)])*1.05,"79", cex = 0.5)

### Hazard Functions : Assuming Exponential Model
library('survival')
exp_model = survreg(Surv(Table_1$survival, !Table_1$cens)~1, dist = 'exponential')

parameters = list('theta' = 1/exp(exp_model$coefficients)); surv_model = 'exponential'

pat_table = matrix(Data[id == j,], ncol = 3)
y = Y[id == j]
cond_dens = Vectorize(h(mean_params, cov_params, parameters, surv_model, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)
norm_probs = probs/sum(probs)
hazard = norm_probs/(1-cumsum(norm_probs))*100

plot(eval_T[eval_T<7],(hazard)[eval_T<7], type = "l", ylim = c(0.1,0.7), xlim = c(5,7), axes = FALSE)
axis(side = 1, at = seq(5,7,0.5), cex.axis = 0.75); axis(side = 2, cex.axis = 0.75)
text(5.2,0.675,"59", cex = 0.5)


y[8] = 69
cond_dens = Vectorize(h(mean_params, cov_params, parameters, surv_model, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)
norm_probs = probs/sum(probs)
hazard = norm_probs/(1-cumsum(norm_probs))*100
lines(eval_T[eval_T<7],hazard[eval_T<7]) 
text(eval_T[1],max(hazard[eval_T<7])*1.05,"69", cex = 0.5)

y[8] = 79
cond_dens = Vectorize(h(mean_params, cov_params, parameters, surv_model, pat_table, y, Covariate))
eval_T = seq(max(pat_table[,1])+0.001,15,0.01)
probs = cond_dens(eval_T)
norm_probs = probs/sum(probs)
hazard = norm_probs/(1-cumsum(norm_probs))*100
lines(eval_T[eval_T<7],hazard[eval_T<7]) 
text(eval_T[1],max(hazard[eval_T<7])*1.05,"79", cex = 0.5)

abline(h = 0.2, lty = 2)
# dev.off()
