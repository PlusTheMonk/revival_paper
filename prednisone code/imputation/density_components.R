g <- function( mean_params, cov_params, pat_table,y,Cov) {
	# Provides a function of the survival time, t,
	# for the likelihood Y | T
	
	g2 <- function(t) {
		k = dim(pat_table)[1]
		X = Cov(pat_table,t)
		S = Sigma(pat_table, t, cov_params)
		mu = X%*%mean_params
		return( (2*pi)^(-k/2)*det(S)^(-1/2)*exp(-t(y - mu)%*% solve(S) %*% (y - mu)/2) )
	}
	
	return(g2)	
}

f <- function(parameters, model = "exponential") { 
	# Assume an Exponential Model with rate parameter theta
	if(model == "exponential") {
    f2 <- function(t) {return(parameters$theta * exp(-parameters$theta * t))}
	}
	if(model == "weibull") {
	  f2 <- function(t) {return(parameters$k/parameters$lambda * (t/parameters$lambda)^(parameters$k-1) * exp(-(t/parameters$lambda)^(parameters$k)) )}
	}
	return(f2)
}

h <- function(mean_params, cov_params, parameters,surv_model, pat_table,y,Cov) {
	# Returns a function that returns value of the joint density 
	# at T = t given the parameter values
	
	h2 <- function(t) {
		t_dens = f(parameters,surv_model)
		y_dens = g(mean_params, cov_params, pat_table,y,Cov) 
		return( t_dens(t) * y_dens(t) )
	}
	
	return(h2)
	
}

#### Imputation  ####

impute_Cov <- function(mean_params, cov_params, parameters, surv_model, pat_table, y, c, Cov) {
	cond_dens = Vectorize(h(mean_params, cov_params, parameters,surv_model, pat_table,y,Cov))

	eval_T = seq(c, c+30, by = 0.01)

	probs = cond_dens(eval_T)

	norm_probs = probs / sum(probs)

	diff = abs(runif(1) - cumsum(norm_probs))
  
  T = eval_T[which(diff == min(diff))]
  
	return(Cov(pat_table,T))
	
}