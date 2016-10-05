## Runs the regression for the application and UI face validity experiment
# 
# Author: mjskay
###############################################################################
library(runjags)
library(tidybayes)
library(metabayes)
library(coda)
library(lme4)
library(plyr)
library(dplyr)
library(ggplot2)
library(rstan)


#------------------------------------------------------------------------------
# THE MODEL
final_model = TRUE

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
    vector<lower=0,upper=1>[n_application] alpha;	#uniform prior		
    ordered[n_acceptable - 1] b0[n_application];
	vector[n_participant] u;
	real<lower=0> tau;
}
transformed parameters {
	real sigma;
	sigma <- 1/sqrt(tau);
}
model {
	#priors for each application
	for (i in 1:n_application) {
        #priors for thresholds
    	b0[i] ~ normal(0, 30);
	}

	#slope of linear model
	b ~ normal(0, 30);

    #power mean
    p ~ normal(0, 2);

	#hyperprior on sd of random effect
	tau ~ gamma(0.001,0.001);

	#priors for random effect: participant
	u ~ normal(0, sigma);

	#core model
	for (i in 1:n) {
    	real mu;
		real ppv_tpr_mean;

		# weighted power mean
		if (p[application[i]] == 0) {
			ppv_tpr_mean <- (ppv[i] ^ alpha[application[i]]) * (tpr[i] ^ (1 - alpha[application[i]]));
		}
		else {
			ppv_tpr_mean <- (alpha[application[i]] * (ppv[i] ^ p[application[i]]) + (1 - alpha[application[i]]) * (tpr[i] ^ p[application[i]])) ^ (1/p[application[i]]);
		}

		#linear model
		mu <- b[application[i]] * ppv_tpr_mean + u[participant[i]];

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

if (!final_model) {
    fit = stan(fit = model, data = data_list, iter = 100, chains = 1, thin = 3)
} else {
    fit = stan(model_code=model_code, data = data_list, iter = 2000, chains = 4)
}

save.image(file=paste0("output/acceptability_ui-model-stan", (if (final_model) "-final" else ""), ".RData"))
#load(file=paste0("output/acceptability_ui-model", (if (final_model) "-final" else ""), ".RData"))

#------------------------------------------------------------------------------
# CHECK FOR CONVERGENCE

print(fit)

pdf(file="output/model-traces.pdf")
traceplot(fit)
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

function ()
foo_data <- ...
rng_seed <- ...
foo <- stan(file = 'foo.stan', data = foo_data, chains = 0)

CL = makeCluster(2, outfile = 'parallel.log')
clusterExport(cl = CL, c("foo_data", "foo", "rng_seed")) 
sflist1 <- parLapply(CL, 1:4, fun = function(cid) {  
          require(rstan)
          stan(fit = foo, data = foo_data, chains = 1, 
                   iter = 2000, seed = rng_seed, chain_id = cid)
    })

fit <- sflist2stanfit(sflist1)
print(fit)
stopCluster(CL)




model = stan(model_code = model_code, chains = 0)
cluster = makeCluster(4, outfile = 'parallel.log')
clusterExport(cluster, c("data_list", "model")) 
fit = parLapply(cluster, 1:4, function(i) {
    require(rstan)
    stan(fit = model, data = data_list, seed = NA,
        iter = 16000, thin=16,
        chains = 1, chain_id = i, 
        refresh = i == 1)
    }) %>%
    sflist2stanfit()


foo <- stan(file = 'foo.stan', data = foo_data, chains = 0)

sflist1 <- parLapply(CL, 1:4, fun = function(cid) {  
          require(rstan)
          stan(fit = foo, data = foo_data, chains = 1, 
                   iter = 2000, seed = rng_seed, chain_id = cid)
    })


