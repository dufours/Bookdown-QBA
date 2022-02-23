# Uploading the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

#I need a R object that will be used in the model text file for the number of observations
a <- dim(cyst)
num_obs <- a[1] 

# Priors for intercept an coefficient of quest ----------------------------

#I need to create R objects that will be used in the model text file to describe the priors
options(scipen=999)           #To avoid use of scientific notation by R
mu_int <- 0.0                 #mu for the intercept    
inv_var_int <- 0.0001         #Inverse variance for the intercept
mu_betaquest <- 0.0           #mu for the coefficient
inv_var_betaquest <- 0.0001   #Inverse variance for the coefficient


# Priors for bias parameters ----------------------------------------------

##################################################################################
##################################################################################
# Below, specify the values that will help describe the prior distributions for
# your bias parameters.
##################################################################################
##################################################################################

#1) I need priors for the log of the unmeasured confounder-disease OR

mu_log_OR_ud <- #TO COMPLETE# 
inv_var_log_OR_ud <- #TO COMPLETE#

#2) I need priors for the prevalence of the confounder in the exposed
library(epiR)

rval1 <- epi.betabuster(#TO COMPLETE# )  
  
rval1$shape1                #View the a shape parameter 
rval1$shape2                ##View the b shape parameter
#plot the prior distribution
curve(dbeta(x, shape1=rval1$shape1, shape2=rval1$shape2), from=0, to=1, 
      main="Prior for prevalence of confounder in exposed", xlab = "Proportion", ylab = "Density")

#3) I need priors for the prevalence of the confounder in the unexposed

rval2 <- epi.betabuster(#TO COMPLETE# )  
  
rval2$shape1                #View the a shape parameter 
rval2$shape2                ##View the b shape parameter 
#plot the prior distribution
curve(dbeta(x, shape1=rval2$shape1, shape2=rval2$shape2), from=0, to=1, 
      main="Prior for prevalence of confounder in unexposed", xlab = "Proportion", ylab = "Density")


# Run Bayesian model using Lash's adjustment for unmeasured confounder ----------------------------

#Generate the .txt file with the model
model_unm_conf_lash <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							test[i] ~ dbern(p[i])
							logit(p[i]) <- int + betaquest*quest[i]  
			}

#===  EXTRA CALCULATION TO GET OR_biased and OR_adj ===#
			OR_biased <- exp(betaquest)
      OR_adj <- OR_biased*( (OR_ud*PU_unexp+1-PU_unexp)/(OR_ud*PU_exp+1-PU_exp) )
			
#===	PRIORS	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		                #FLAT PRIOR
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#FLAT PRIOR
			
			log_OR_ud ~ dnorm(", mu_log_OR_ud,"," , inv_var_log_OR_ud,")   #PRIOR FOR log_OR_ud		
			OR_ud <- exp(log_OR_ud)	                                       #LINK OR_ud TO THE LOG_OR_UD
			
			PU_exp ~ dbeta(", rval1$shape1,",", rval1$shape2, ")		       #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE EXPOSED INDIVIDUALS 
			PU_unexp ~ dbeta(", rval2$shape1,",", rval2$shape2, ")		     #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE UNEXPOSED INDIVIDUALS 	

		}")

#write to temporary text file
write.table(model_unm_conf_lash, file="model_unm_conf_lash.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

##################################################################################
##################################################################################
# Now you need a number of sets of initial values corresponding to the number of 
# Markov chains that will be ran. Use three chains for now.
# Provide values for the three chains. 
##################################################################################
##################################################################################

#Now I need a number of sets of initial values corresponding to the number of Markov chains that will be ran. I plan to use 3 chains, thus 3 sets of initial values
#Initializing values for 3 chains
inits <- list(
  
  list(#TO COMPLETE# ),            
  list(#TO COMPLETE# ),
  list(#TO COMPLETE# )
    
)

# We have everything, we can now run the model using R2OpenBUGS
library(R2OpenBUGS)

#########################
#########################
#Set number of iterations
#########################
#########################

#Set number of iterations

niterations <- #TO COMPLETE#
  
#Run the Bayesian model
bug.out.unm_conf_lash <- bugs(data=cyst, 
                      inits=inits, 
                      parameters.to.save=c("int", "betaquest", "OR_biased", "OR_adj"),  
                      n.iter=niterations, 
                      n.burnin=0, 
                      n.thin=1,
                      n.chains=3, 
                      model.file="model_unm_conf_lash.txt", 
                      debug=TRUE,
                      DIC=FALSE)






# MCMC diagnostic ---------------------------------------------------------

#####################################################
#####################################################
# Used the bug.out.unm_conf_lash R object to generate 
# the diagnostic plots.
#####################################################
#####################################################

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)

#TO COMPLETE#


# Applying burn in ---------------------------------

#####################################################
#####################################################
# Set up a burn in period, 
# discard values generated by the MCMC during this period, 
# and combine the post burn in values of the Markov chains.
#####################################################
#####################################################

# Now we need to set up a burn in period and discard values generated by the MCMC during this period
# Set burn-in to 1000

burnin <- #TO COMPLETE#
  
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_unm_conf_lash <- data.frame(rbind(bug.out.unm_conf_lash$sims.array[k,1,],bug.out.unm_conf_lash$sims.array[k,2,],bug.out.unm_conf_lash$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_unm_conf_lash <- round(
  t(
    apply(X=estimates_unm_conf_lash, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_unm_conf_lash

# Plot the prior and posterior distributions for betaquest
std <- sqrt(1/inv_var_betaquest)
plot(density(x=estimates_unm_conf_lash$betaquest), 
     main="Quest coefficient",
     xlab="Value", ylab="Density",
)
curve(dnorm(x, 
            mean=(mu_betaquest), 
            sd=std
),
lty=2,
col="red",
add=TRUE
)