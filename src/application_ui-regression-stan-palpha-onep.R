## Runs the regression for the application and UI face validity experiment
# 
# Author: mjskay
###############################################################################
library(tidybayes)
library(plyr)
library(dplyr)
library(rstan)
library(ggmcmc)
library(parallel)

memory.limit(8096)

#------------------------------------------------------------------------------
# THE MODEL
final_model = TRUE

model = metastan(
    data = {
        n: int(lower=1)
        n_application: int(lower=1)
        application[n]: int(lower=1,upper=n_application)
        n_acceptable: int(lower=2)
        acceptable[n]: int(lower=1,upper=n_acceptable)
        n_participant: int(lower=1)
        participant[n]: int(lower=1,upper=n_participant)
        ppv[n]: real
        tpr[n]: real
    },
    parameters = {
        b: vector[n_application]
        p: vector[n_application]
        alpha_mu: vector(lower=0,upper=1)[n_application]
        alpha_v: vector(lower=0)[n_application]
        u_alpha[n_participant]: vector(lower=0,upper=1)[n_application]
        b0[n_application]: ordered[n_acceptable - 1]
        u_b0: vector[n_participant]
        u_b: vector[n_participant]
        sigma_b: real(lower=0)
        sigma_b0: real(lower=0)
    },
    model = {
    	#priors for each application
    	for (i in 1:n_application) {
            #priors for thresholds
        	b0[i] ~ normal(0, 10);
    	}
        
    	#slope of linear model
    	b ~ normal(0, 10);
        
        #power mean
        p ~ normal(0, 1);
        
    	#hyperprior on sd of random effects
    	sigma_b0 ~ gamma(1, 1);
    	sigma_b ~ gamma(1, 1);
        
    	#priors for participant effect: intercept
    	u_b0 ~ normal(0, sigma_b0);
        
    	#priors for participant effect: slope
    	u_b ~ normal(0, sigma_b);
        
    	#priors for participant effect: individual alpha
    	alpha_v ~ gamma(1, 1);
    	for (k in 1:n_participant) {
    		u_alpha[k] ~ beta(alpha_mu %*% alpha_v, (1 - alpha_mu) %*% alpha_v);
    	}
        
    	#core model
    	for (i in 1:n) {
        	mu: real
    		ppv_tpr_mean: real
            
    		# weighted power mean
    		if (p[application[i]] == 0) {
    			ppv_tpr_mean <- (ppv[i] ^ u_alpha[participant[i],application[i]]) * (tpr[i] ^ (1 - u_alpha[participant[i],application[i]]));
    		}
    		else {
    			ppv_tpr_mean <- (u_alpha[participant[i],application[i]] * (ppv[i] ^ p[application[i]]) + (1 - u_alpha[participant[i],application[i]]) * (tpr[i] ^ p[application[i]])) ^ (1/p[application[i]]);
    		}
            
    		#linear model
    		mu <- (b[application[i]] + u_b[participant[i]]) * ppv_tpr_mean + u_b0[participant[i]];
            
    		#response
    	    acceptable[i] ~ ordered_logistic(mu, b0[application[i]]);
      	}
})

model_code = "
data {
    int<lower=1> n;
    int<lower=1> n_application;
    int<lower=1,upper=n_application> application[n];
    int<lower=2> n_acceptable;
    int<lower=1,upper=n_acceptable> acceptable[n];
    int<lower=1> n_participant;
    int<lower=1,upper=n_participant> participant[n];
	real ppv[n];
	real tpr[n];
}
parameters {
    vector[n_application] b;
    vector[n_application] p;
    vector<lower=0,upper=1>[n_application] alpha_mu;		
    vector<lower=0>[n_application] alpha_v;
    vector<lower=0,upper=1>[n_application] u_alpha[n_participant];
    ordered[n_acceptable - 1] b0[n_application];
	vector[n_participant] u_b0;
	vector[n_participant] u_b;
	real<lower=0> sigma_b;
	real<lower=0> sigma_b0;
}
model {
	#priors for each application
	for (i in 1:n_application) {
        #priors for thresholds
    	b0[i] ~ normal(0, 10);
	}

	#slope of linear model
	b ~ normal(0, 10);

    #power mean
    p ~ normal(0, 1);

	#hyperprior on sd of random effects
	sigma_b0 ~ gamma(1, 1);
	sigma_b ~ gamma(1, 1);

	#priors for participant effect: intercept
	u_b0 ~ normal(0, sigma_b0);

	#priors for participant effect: slope
	u_b ~ normal(0, sigma_b);

	#priors for participant effect: individual alpha
	alpha_v ~ gamma(1, 1);
	for (k in 1:n_participant) {
		u_alpha[k] ~ beta(alpha_mu .* alpha_v, (1 - alpha_mu) .* alpha_v);
	}

	#core model
	for (i in 1:n) {
    	real mu;
		real ppv_tpr_mean;

		# weighted power mean
		if (p[application[i]] == 0) {
			ppv_tpr_mean <- (ppv[i] ^ u_alpha[participant[i],application[i]]) * (tpr[i] ^ (1 - u_alpha[participant[i],application[i]]));
		}
		else {
			ppv_tpr_mean <- (u_alpha[participant[i],application[i]] * (ppv[i] ^ p[application[i]]) + (1 - u_alpha[participant[i],application[i]]) * (tpr[i] ^ p[application[i]])) ^ (1/p[application[i]]);
		}

		#linear model
		mu <- (b[application[i]] + u_b[participant[i]]) * ppv_tpr_mean + u_b0[participant[i]];

		#response
	    acceptable[i] ~ ordered_logistic(mu, b0[application[i]]);
  	}
}
"

#------------------------------------------------------------------------------
# THE DATA.
source("src/application_ui-load_data.R")

#build data list
data_list = df %>%
    filter(!is.na(acceptable)) %>%	#drop NA responses (there are only 2 anyway)
    select(acceptable, application, participant, tpr, ppv) %>%
    compose_data()

#------------------------------------------------------------------------------
# FIT THE MODEL

model = stan(model_code = model_code, data = data_list, chains = 0)

if (!final_model) {
    iterations = 1000
    thin = 1
} else {
    iterations = 40000
    thin = 4
}

cluster = makeCluster(4, outfile = 'parallel.log')
clusterExport(cluster, c("data_list", "model", "iterations", "thin")) 
fit = parLapply(cluster, 1:4, function(i) {
                require(rstan)
                stan(fit = model, data = data_list, seed = NA,
                        iter = iterations, thin = thin,
                        chains = 1, chain_id = i, 
                        refresh = i == 1)
            }) %>%
    sflist2stanfit()


save.image(file=paste0("output/acceptability_ui-model-stan", (if (final_model) "-final" else ""), ".RData"))
load(file=paste0("output/acceptability_ui-model-stan", (if (final_model) "-final" else ""), ".RData"))

#------------------------------------------------------------------------------
# CHECK FOR CONVERGENCE

print(fit)

pdf(file="output/model-traces.pdf")
rstan:::traceplot(fit)
dev.off()

ggmcmc(ggs(fit), file="output/model-autocorr.pdf", plot="ggs_autocorrelation")

ggmcmc(ggs(fit), file="output/model-density.pdf", plot="ggs_density")

#------------------------------------------------------------------------------
# EXTRACT THE RESULTS

# Convert coda-object coda_samples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmc_chain[ stepIdx , paramIdx ]
fit = apply_prototypes(fit, df)

mcmc_chain = as.matrix(coda_samples) %>%
    apply_prototypes(df)
#drop other (often large) fit objects; no longer needed
rm(jags_fit, coda_samples)

#subset to preferred model

#save
#save.image(file=paste0("output/acceptability_ui-model-small", (if (final_model) "-final" else ""), ".RData"))




library(parallel)
library(rstan)




cluster = makeCluster(4, outfile = 'parallel.log')
clusterExport(cluster, c("data_list", "model")) 
fit = parLapply(cluster, 1:4, function(i) {
    require(rstan)
    stan(fit = model, data = data_list, seed = NA,
        iter = 40000, thin=4,
        chains = 1, chain_id = i, 
        refresh = i == 1)
    }) %>%
    sflist2stanfit()
