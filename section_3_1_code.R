lambda= 1/10

t = c(0,1,2,3)
y = c(6.0,4.5,5.4,4.0)

S = diag(4)+exp(-abs(outer(t,t,"-")))

mu <- function(beta) {
  mu_T <- function(T) {
    s = T-t
    return(beta*s/(1+s))
  }
  return(mu_T)
}

res <- function(beta) {
  res_T <- function(T) {
    e = y-mu(beta)(T)
    return(exp(-t(e)%*%solve(S,e)/2)*dexp(T,rate = lambda))
  }
}

eval_T = seq(3.01,30,0.01)

### Left plot
# png("/Users/walterdempsey/Documents/stat/research/joint_models/censoring/survival_models/conditional_density.png", width = 6.5,height = 3, units = "in", res = 300)
op <- par(mfrow = c(1,2),
          oma = c(2,1,0,0) + 0.1,
          mar = c(0,1,1,1) + 0.1)

probs = Vectorize(res(0))(eval_T)
norm_probs = probs/sum(probs)*100

plot(eval_T, (norm_probs), type = "l", ylim = c(0.0,0.1), xlim = c(3,30), axes = FALSE, ylab = "", xlab = "")
axis(side = 1, at = seq(5,30,5), cex.axis = 0.75);axis(side = 2, at = seq(0,0.1,0.05), cex.axis = 0.75)

beta_seq = seq(1/3,2, length.out = 6)

for(beta in beta_seq) {
  probs = Vectorize(res(beta))(eval_T)
  norm_probs = probs/sum(probs)*100
  lines(eval_T, (norm_probs), type = "l")
}

text(7,0.01,expression(paste(beta, "= 2.0")), cex = 0.75)
text(7,0.1,expression(paste(beta, "= 0.0")), cex = 0.75)

### Right plot
probs = Vectorize(res(4))(eval_T)
norm_probs = probs/sum(probs)*100

plot(eval_T, (norm_probs-norm_probs[1]), type = "l", ylim = c(0.0,0.5), xlim = c(3,30), axes = FALSE, ylab = "", xlab = "")
axis(side = 1, at = seq(5,30,5), cex.axis = 0.75);axis(side = 2, at = seq(0,0.5,0.1), cex.axis = 0.75)

beta_seq = seq(5,8)

for(beta in beta_seq) {
  probs = Vectorize(res(beta))(eval_T)
  norm_probs = probs/sum(probs)*100
  lines(eval_T, (norm_probs-norm_probs[1]), type = "l")
}

text(7,0.5,expression(paste(beta, "= 8.0")), cex = 0.75)
text(27.5,0.055,expression(paste(beta, "= 4.0")), cex = 0.75)

# dev.off()