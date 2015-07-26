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


#------------------------------------------------------------------------------
# THE MODEL
final_model = FALSE

model = metajags({
	#core model
	for (i in 1:n) {
		# linear model
		mu[i] <- 
			#weighted power mean M_p
			ifelse(p[application[i]] == 0,
                b[application[i]] * ((ppv[i] ^ alpha[application[i]]) * (tpr[i] ^ (1 - alpha[application[i]]))),
                b[application[i]] * (alpha[application[i]] * (ppv[i] ^ p[application[i]]) + (1 - alpha[application[i]]) * (tpr[i] ^ p[application[i]])) ^ one_over_p[application[i]]
            ) + u[participant[i]]   #participant effect
        
		# logits for thresholds between adjacent levels
        logit(Q[i,1]) <- b0[application[i],1] - mu[i] 
        P[i,1] <- Q[i,1]
        for (j in 2:(n_acceptable - 1)) {
            logit(Q[i,j]) <- b0[application[i],j] - mu[i]
            P[i,j] <- Q[i,j] - Q[i,j-1]
        }
        P[i,n_acceptable] <- 1 - Q[i,n_acceptable - 1] 

        #observation
        acceptable[i] ~ dcat(P[i,])
  	}

	#priors for each scenario
  	for (i in 1:n_application) {
        #priors for thresholds
        for (j in 1:(n_acceptable - 1)) {
	            b0_unsorted[i,j] ~ dnorm(0, 1.0E-3)
        }
        b0[i,1:(n_acceptable - 1)] <- sort(b0_unsorted[i,1:(n_acceptable - 1)])

		#slope of linear model
    	b[i] ~ dnorm(0, 1.0E-3)

		#weight on ppv for mean of ppv and tpr
		alpha[i] ~ dunif(0, 1)

        #power mean
        p[i] ~ dnorm(0, 0.25)
        #hack to prevent divide by zero error in core model when p == 0
        #(Note the core model doesn't actually use the value of one_over_p
        #when p == 0, so we can set it to anything here). 
        one_over_p[i] <- 1/ifelse(p[i] == 0, 1, p[i])
  	}

	#priors for random effect: participant
	for (k in 1:n_participant) {
		u[k] ~ dnorm(0, tau)
	}

	#hyperprior on sd of random effect
	tau ~ dgamma(0.001,0.001) 
})


#------------------------------------------------------------------------------
# THE DATA.
source("src/application_ui-load_data.R")

#build data list
data_list = df %>%
    select(acceptable, application, participant, tpr, ppv) %>%
    compose_data()

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

#get maximum likelihood estimates of some model parameters from the p=1 (arithmetic mean)
#version of the model (since it is easiest to calculate) using ordinal
#logistic regression. We will use these as starting points for the model.
models = dlply(df, ~ application, function(df) {
        clmm(acceptable ~ ppv + tpr + (1|participant), data=df)
    })
inits_list = ldply(models, function(m) {
        b.ppv = coef(m)["ppv"]
        b.tpr = coef(m)["tpr"]
        data.frame(
            b = b.ppv + b.tpr,
            alpha = b.ppv / (b.ppv + b.tpr),
            p = 1   #inits are for power mean with p=1 (arithmetic mean)
        )
    }) %>%
    select(-application) %>%
    as.list()

#inits_list=list()
inits_list$b0_unsorted = laply(models, function(m) m$alpha)



#------------------------------------------------------------------------------
# FIT THE MODEL

inits_function = function(chain) {
    c(inits_list, .RNG.name = c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")[[chain]])
}

parameters = c("b0" , "b", "alpha", "p")  # The parameter(s) to be monitored.
if (!final_model) {
    jags_fit = run.jags(code(model), data=data_list, monitor=parameters, inits=list(inits_list, inits_list), modules="glm",
        method="parallel")
} else {
#    jags_fit = autorun.jags(code(model), data=data_list, monitor=parameters, initlist=inits_list,
#        method="parallel", thin.sample=TRUE)    
    jags_fit = run.jags(code(model), data=data_list, monitor=parameters, initlist=inits_list,
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
mcmc_chain = as.matrix(coda_samples) %>%
    apply_prototypes(df)
#drop other (often large) fit objects; no longer needed
rm(jags_fit, coda_samples)

#subset to preferred model

#save
#save.image(file=paste0("output/acceptability_ui-model-small", (if (final_model) "-final" else ""), ".RData"))

