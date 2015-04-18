## Runs the regression for the application and UI face validity experiment
# 
# Author: mjskay
###############################################################################
library(runjags)
library(tidybayes)
library(metajags)
library(coda)
library(lme4)
library(plyr)
library(dplyr)
library(ggplot2)


#------------------------------------------------------------------------------
# THE MODEL
final_model = TRUE

model = metajags_model({
	#core model
	for (i in 1:n) {
    	acceptable_b[i] ~ dbern(mu[i])

		# linear model with logit link
		logit(mu[i]) <- b0[application[i]] +
			#harmonic mean (M_-1, weighted F-measure)
			equals(model_number, 1) * (b[application[i]]/(alpha[application[i]]/ppv[i] + (1-alpha[application[i]])/tpr[i])) +
			#geometric mean (M_0, weighted G-measure)
			equals(model_number, 2) * (b[application[i]]*(pow(ppv[i], alpha[application[i]]) * pow(tpr[i],(1-alpha[application[i]])))) +
			#arithmetic mean (M_1) 
			equals(model_number, 3) * (b[application[i]]*(alpha[application[i]]*ppv[i] + (1-alpha[application[i]])*tpr[i])) +
			u[participant[i]]
  	}

	#priors for each scenario
  	for (i in 1:n_application) {
		#intercept
	  	b0[i] ~ dnorm(0, 1.0E-3)

		#slope of linear model
    	b[i] ~ dnorm(0, 1.0E-3)

		#weight on ppv for mean of ppv and tpr
		alpha[i] ~ dunif(0, 1)
  	}

	#priors for random effect: participant
	for (k in 1:n_participant) {
		u[k] ~ dnorm(0, tau)
	}

	#hyperprior on sd of random effect
	tau ~ dgamma(0.001,0.001)

	#model comparison hyperprior
	model_prob[1] <- 1/3
	model_prob[2] <- 1/3
	model_prob[3] <- 1/3
	model_number ~ dcat(model_prob[])
})


#------------------------------------------------------------------------------
# THE DATA.
source("src/application_ui-load_data.R")

#build data list
data_list = df %>%
    select(acceptable_b, application, participant, tpr, ppv) %>%
    compose_data()

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

#get maximum likelihood estimates of model parameters from the p=1 (arithmetic mean)
#version of the model (since it is easiest to calculate) using traditional mixed-effects
#logistic regression. We will use these as starting points for the model.
inits_list = ddply(df, ~ application, function(df) { 
        m = glmer(acceptable_b ~ ppv + tpr + (1|participant), data=df, family=binomial)
        b.ppv = fixef(m)["ppv"]
        b.tpr = fixef(m)["tpr"]
        data.frame(
	            b0 = fixef(m)["(Intercept)"],
	            b = b.ppv + b.tpr,
	            alpha = b.ppv / (b.ppv + b.tpr)
        )
    }) %>%
    select(-application)


#------------------------------------------------------------------------------
# FIT THE MODEL

parameters = c("b0" , "b", "alpha", "model_number")  # The parameter(s) to be monitored.
if (!final_model) {
    jags_fit = run.jags(model$code, data=data_list, monitor=parameters, initlist=inits_list, 
        method="parallel")
} else {
#    jags_fit = autorun.jags(model$code, data=data_list, monitor=parameters, initlist=inits_list,
#        method="parallel", thin.sample=TRUE)    
    jags_fit = run.jags(model$code, data=data_list, monitor=parameters, initlist=inits_list,
	    adapt=1000, burnin = 50000, sample = 10000, thin=15,
        method="parallel")
}

#save.image(file=paste0("output/acceptability_ui-model", (if (final_model) "-final" else ""), ".RData"))
#load(file=paste0("output/acceptability_ui-model", (if (final_model) "-final" else ""), ".RData"))

#------------------------------------------------------------------------------
# CHECK FOR CONVERGENCE

plot(jags_fit, file="output/model-params.pdf")
summary(jags_fit)

pdf(file="output/model-autocorr.pdf")
coda_samples = as.mcmc.list(jags_fit)
autocorr.plot(coda_samples, ask=FALSE)
dev.off()

#------------------------------------------------------------------------------
# EXTRACT THE RESULTS

# Convert coda-object coda_samples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmc_chain[ stepIdx , paramIdx ]
mcmc_chain = as.matrix(coda_samples)
#drop other (often large) fit objects; no longer needed
rm(jags_fit, coda_samples)

#subset to preferred model
model_number = mcmc_chain[,"model_number"]
p_model = prop.table(table(model_number))
best_model = which(p_model == max(p_model))
best_model_chain = mcmc_chain[model_number == best_model,]  %>%
    apply_prototypes(df)

#save
save.image(file=paste0("output/acceptability_ui-model-small", (if (final_model) "-final" else ""), ".RData"))

