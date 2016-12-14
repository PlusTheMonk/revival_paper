simulate_survival <- function(num_sims = 1, parameters, X, censoring_time = 0, baseline = 'exponential') {
	## Simulate Survival Times Given a Proportional Hazards Model with a Chosen Baseline Hazard Rate
	## When the Covariate, X, are time-independent
	
	if(baseline == 'exponential') {
		exp_rate = parameters$lambda*exp(-parameters$beta*X)
		return(censoring_time+rexp(num_sims, rate = exp_rate))
	} else if(baseline == 'weibull') {
		alpha = exp(-parameters$beta*X)
		lambda = parameters$lambda
		k = parameters$k
		lower_bound = 1-exp(-(censoring_time/k)^lambda*(1/alpha))
		return(qweibull(1-(1-runif(num_sims, min = lower_bound))^(alpha), shape = lambda, scale = k))
	} else if(baseline == 'piecewise') {
		breaks = parameters$breaks[parameters$breaks >= censoring_time]
		values = parameters$values[parameters$breaks >= censoring_time]
		alpha = exp(-parameters$beta*X)
		
		T = censoring_time + rexp(num_sims, rate = values[1])			
		T[T > breaks[1]] = breaks[1]
		if(length(breaks) > 1) {
			for (i in 2:length(breaks)) {
					T = T + rexp(num_sims, rate = values[i]*alpha)*(T >= breaks[i-1])
					T[T> breaks[i]] = breaks[i]					
			}			
		}
		return(T)				
	} else if(baseline == 'prop_haz') {		
		times = parameters$times[parameters$times > censoring_time]
		deaths = parameters$deaths[parameters$times > censoring_time]
		T = rep(max(parameters$times[parameters$deaths]),num_sims) # initialize to final death time
		
		deaths = deaths[order(times)]; times = times[order(times)]

		if (censoring_time > max(parameters$times[parameters$deaths])) {
			T = rep(censoring_time,num_sims)
		} else {		
			for (k in 1:num_sims) {
				for (t in times[deaths]) {
					risk = which(times > t)
					death = which((times == t) & (deaths))
					hazard = sum(exp(-as.matrix(X[death,])%*%beta))/sum(exp(-as.matrix(X[risk,])%*%beta))
					if(runif(1) < hazard) {T[k] = t; break}
				}	
			}
		}		
		return(T)
	} else {
		cat("Error: Incorrect Specification of Baseline Model /n")
	}
	
	
}