library('VGAM')

draw_fit_gg = function(converted, total, prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I, passage) {
		ratio = converted/total

    # Compute likelihoods
	# Background likelihood (assuming discrete uniform distribution)
    likelihood_H = 1 / total
    # Alternative background likelihood, bounded
    #likelihood_H = (ratio <= mu_U & ratio >= mu_M) / ((mu_U - mu_M) * total)
    # Second alternative background likelihood, with shoulders
    #likelihood_H = ((ratio/mu_M) * (ratio < mu_M) + (ratio <= mu_U & ratio >= mu_M) + (1 - ratio)/(1 - mu_U) * (ratio > mu_U)) / ((mu_U/2 - mu_M/2 + 1/2) * total)
	# Third alternative background
	#likelihood_H = dnormp(converted, p=8, mu=500, sigmap = 300)       # mu_H =; sigmap =
	# Fourth alternative background
	#tot<-(ratio/mu_M) * (ratio < mu_M) + (ratio <= mu_U & ratio >= mu_M)  + (1 - ratio)/(1 - mu_U) * (ratio > mu_U)
	#tot[which(!(ratio <= mu_U | ratio >= mu_M))]<-tot[which(!(ratio <= mu_U | ratio >= mu_M))]-1
	#likelihood_H<-tot/((mu_U/2 - mu_M/2 + 1/2) * total)



	likelihood_U = dbetabinom(converted, total, mu_U, rho_U)
    likelihood_M = dbetabinom(converted, total, mu_M, rho_M)
    likelihood_I = dbetabinom(converted, total, mu_I, rho_I)
    L = sum(log(likelihood_H + likelihood_U + likelihood_M + likelihood_I))

    
    # Compute expected_likelihoods
    expected_likelihood_H = prior_H * likelihood_H
    expected_likelihood_U = prior_U * likelihood_U
    expected_likelihood_M = prior_M * likelihood_M
    expected_likelihood_I = prior_I * likelihood_I
    expected_likelihood_total = expected_likelihood_H + expected_likelihood_U + expected_likelihood_M + expected_likelihood_I
				
		plt<-data.frame(
		        beta = sort(ratio), 
		        unmethylated = expected_likelihood_U[order(ratio)] * total[order(ratio)], 
                background = expected_likelihood_H[order(ratio)] * total[order(ratio)],
                methylated = expected_likelihood_M[order(ratio)] * total[order(ratio)],
                intermediate = expected_likelihood_I[order(ratio)] * total[order(ratio)],
                total = expected_likelihood_total[order(ratio)]* total[order(ratio)])

	ggplot()+
		  geom_histogram(aes(x=beta, y=..density..),plt,bins=100, color="grey85",fill="grey85")+
		  geom_line(aes(beta,background), plt, color="grey40", size=0.75)+
		  geom_line(aes(beta,unmethylated), plt, color="green", size=0.75)+
		  geom_line(aes(beta,methylated), plt, color="blue", size=0.75)+
		  geom_line(aes(beta,intermediate), plt, color="red", size=0.75)+
		  geom_line(aes(beta,total), plt, color="black", size=0.75)+
		  annotate("text", x=0.5, y=3, label=paste(passage, "\n Log Likelihood: ", round(L,2), sep=""), color="grey30", size=5)+theme_bw()+
		  theme(axis.text = element_text(size =10),
		        axis.title = element_text(size =12))+ xlab("DNAm Beta")+ylim(0,5)

}

draw_fit_params_gg = function(converted, total, model, passage) {
	draw_fit_gg(converted, total, model$prior_H, model$prior_U, model$prior_M, model$prior_I, model$mu_U, model$rho_U, model$mu_M, model$rho_M, model$mu_I, model$rho_I, passage)
}



# Maximising log-likelihood as a function of beta-binomial params for one model (see http://varianceexplained.org/r/mixture-models-baseball/)
maximise_betabinom_param = function(converted, total, membership_probability, mu, rho) {
  ll = function(mu, rho) {
    -sum(membership_probability * dbetabinom(converted, total, mu, rho, log = TRUE))
  }
  stats4::coef(stats4::mle(ll, start = list(mu = mu, rho = rho), method = "L-BFGS-B", lower=c(0.001,0.001),upper=c(0.999,0.999)))
}


em = function(converted, total, prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I) {	
	last_L = 0
	L = 10 #1e-10
	ratio = converted/total
	print(paste(c("loglik","Prior(H)","Prior(U)","Prior(M)","Prior(I)", "Mu(U)", "Rho(U)", "MU(M)", "Rho(M)", "mu_I", "rho_I")))
	print(paste(c("XXX", prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I), sep='\t'))
			
	# Convergence determined on expected_likelihood change
	while (abs(L - last_L) >= 10 ) { #(abs(0.001 * last_L))
		last_L = L
		
		# Compute likelihoods
		# Background likelihood (assuming discrete uniform distribution)
		likelihood_H = 1 / total
		# Alternative background likelihood, bounded
		#likelihood_H = (ratio <= mu_U & ratio >= mu_M) / ((mu_U - mu_M) * total)
		# Second alternative background likelihood, with shoulders
		#likelihood_H = ((ratio/mu_M) * (ratio < mu_M) + (ratio <= mu_U & ratio >= mu_M) + (1 - ratio)/(1 - mu_U) * (ratio > mu_U)) / ((mu_U/2 - mu_M/2 + 1/2) * total)
		# Third alternative background
		#likelihood_H = dnormp(converted, p=8, mu=500, sigmap = 300)       # mu_H =; sigmap =
		# Fourth alternative background
		#tot<-(ratio/mu_M) * (ratio < mu_M) + (ratio <= mu_U & ratio >= mu_M)  + (1 - ratio)/(1 - mu_U) * (ratio > mu_U)
		#tot[which(!(ratio <= mu_U | ratio >= mu_M))]<-tot[which(!(ratio <= mu_U | ratio >= mu_M))]-1
		#likelihood_H<-tot/((mu_U/2 - mu_M/2 + 1/2) * total)



		likelihood_U = dbetabinom(converted, total, mu_U, rho_U)
		likelihood_M = dbetabinom(converted, total, mu_M, rho_M)
		likelihood_I = dbetabinom(converted, total, mu_I, rho_I)
		
		# Compute expected likelihoods 
		expected_likelihood_H = prior_H * likelihood_H
		expected_likelihood_U = prior_U * likelihood_U
		expected_likelihood_M = prior_M * likelihood_M
		expected_likelihood_I = prior_I * likelihood_I
		expected_likelihood = expected_likelihood_H + expected_likelihood_M + expected_likelihood_U + expected_likelihood_I

		# Compute membership probabilities (E step)
    	L = sum(log(likelihood_H + likelihood_U + likelihood_M + likelihood_I))
		membership_probability_H = expected_likelihood_H / expected_likelihood
		membership_probability_U = expected_likelihood_U / expected_likelihood
		membership_probability_M = expected_likelihood_M / expected_likelihood
		membership_probability_I = expected_likelihood_I / expected_likelihood
		
		# M step (see https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Gaussian_mixture)
		# Re-evaluate priors
		prior_H = mean(membership_probability_H)
		prior_U = mean(membership_probability_U)
		prior_M = mean(membership_probability_M)
		prior_I = mean(membership_probability_I)

		# Optimize BetaBinomial params 
		mu_rho = maximise_betabinom_param(converted, total, membership_probability_U, mu_U, rho_U)
		mu_U = mu_rho[[1]]
		rho_U = mu_rho[[2]]
		
		mu_rho = maximise_betabinom_param(converted, total, membership_probability_M, mu_M, rho_M)
		mu_M = mu_rho[[1]]
		rho_M = mu_rho[[2]]
		
		mu_rho = maximise_betabinom_param(converted, total, membership_probability_I, mu_I, rho_I)
		mu_I = mu_rho[[1]]
		rho_I = mu_rho[[2]]
		
		# Print results
		print(round(c(L, prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I),2))
	}
	
	list(prior_H = prior_H, prior_U = prior_U, prior_M = prior_M, prior_I = prior_I, mu_U = mu_U, rho_U = rho_U, mu_M = mu_M, rho_M = rho_M, mu_I = mu_I, rho_I = rho_I)
}




test_em = function(prior_H=0.25, prior_U=0.25, prior_M=0.25, prior_I=0.25, mu_U = .90, rho_U = .01, mu_M = .1, rho_M = .01, mu_I = 0.6, rho_I = .01) {
	# Generate dataset
	total = 1000
	min_count = 100
	
	hemimethylated_counts = round(rexp(prior_H * total, rate=.01)) + min_count
	hemimethylated_converted = round(runif(length(hemimethylated_counts)) * hemimethylated_counts)
	#Alternative background distribution, bounded
	# hemimethylated_converted = round(runif(length(hemimethylated_counts), min = mu_M, max=mu_U) * hemimethylated_counts)
	methylated_counts = round(rexp(prior_M * total, rate=.01)) + min_count
	methylated_converted = rbetabinom(length(methylated_counts), methylated_counts, mu_M, rho_M)
	unmethylated_counts = round(rexp(prior_U * total, rate=.01)) + min_count
	unmethylated_converted = rbetabinom(length(unmethylated_counts), unmethylated_counts, mu_U, rho_U)
	intermediate_counts = round(rexp(prior_I * total, rate=.01)) + min_count
	intermediate_converted = rbetabinom(length(intermediate_counts), intermediate_counts, mu_I, rho_I)
	counts = c(hemimethylated_counts, methylated_counts, unmethylated_counts, intermediate_counts)
	converted = c(hemimethylated_converted, methylated_converted, unmethylated_converted, intermediate_converted)
	
	
	# Run test
	res = em(converted, counts, prior_H, prior_U, prior_M, prior_I, .1, .1, .90, .03, .5, .05)
	print(res)
	draw_fit_params_gg(converted, counts, res, 42)
}



# Comparing likelihoods
fit_params_L = function(converted, total, model) {
  fit_L(converted, total, model$prior_H, model$prior_U, model$prior_M, model$prior_I, model$mu_U, model$rho_U, model$mu_M, model$rho_M, model$mu_I, model$rho_I)
}

fit_L = function(converted, total, prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I) {
  ratio = converted/total
  likelihood_H = 1 / total
  likelihood_U = dbetabinom(converted, total, mu_U, rho_U)
  likelihood_M = dbetabinom(converted, total, mu_M, rho_M)
  likelihood_I = dbetabinom(converted, total, mu_I, rho_I)
    
  L = sum(log(likelihood_H + likelihood_U + likelihood_M + likelihood_I))
  L
}




passage_threshold_params = function(converted, total, model) {
	passage_threshold(converted, total, model$prior_H, model$prior_U, model$prior_M, model$prior_I, model$mu_U, model$rho_U, model$mu_M, model$rho_M, model$mu_I, model$rho_I)
}

passage_threshold = function(converted, total, prior_H, prior_U, prior_M, prior_I, mu_U, rho_U, mu_M, rho_M, mu_I, rho_I) {
		# Compute likelihoods
		likelihood_H = prior_H * 1/total
		likelihood_U = prior_U * dbetabinom(converted, total, mu_U, rho_U)
		likelihood_M = prior_M * dbetabinom(converted, total, mu_M, rho_M)
		likelihood_I = prior_I * dbetabinom(converted, total, mu_I, rho_I)
		likelihood_total = likelihood_H + likelihood_M + likelihood_U + likelihood_I

		ratio = converted/total
		max(likelihood_I)/max(likelihood_H)
		}