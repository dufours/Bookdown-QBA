# Uploading and exploring the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

#Create a "a" object with number of rows and columns
a <- dim(cyst)
#Save and check the first value of the "a" object
num_obs <- a[1]

# Priors for intercept and coefficient -------------------------------------------------------------
options(scipen=999) # To avoid use of scientific notation by R

mu_int <- 0.0         #mu for the intercept    
inv_var_int <- 0.0001 #Inverse variance for the intercept
mu_betaquest <- 0.0   #mu for the coefficient
inv_var_betaquest <- 0.0001 #Inverse variance for the coefficient


# Values for the bias parameters ------------------------------------------
##########################################################################################################
# Rather than proposing prior distributions, you could first use exact values (deterministic approach) ###
##########################################################################################################

###TO COMPLETE###

##########################################################################################
### Specify the measurement error part of the model's likelihood function              ###
### Add the values for Se and Sp                                                       ###
##########################################################################################

model_exp_non_dif1 <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							
							#RESPONSE PART
							test[i] ~ dbern(P[i])
							logit(P[i]) <- int + betaquest*expo_true[i] 
							
							#MEASUREMENT ERROR PART
							
							
							### TO COMPLETE ###
					
			}

#===   EXTRA CALCULATION TO GET OR ===#
			ORquest <- exp(betaquest)
			
#===	PRIOR	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		#PRIOR INTERCEPT
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#PRIOR COEFFICIENT
			Pexpo ~ dbeta(1.0, 1.0)                                     #FLAT PRIOR ON PREVALENCE OF TRUE EXPOSURE
			
			
			### TO COMPLETE ###
			
			}")

#write to temporary text file
write.table(model_exp_non_dif1, file="model_exp_non_dif1.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

# Generate initial values -------------------------------------------------

# In this case, I do not need initial values for the bias parameters, since we specified exact values for them (vs distributions)
# But I do need a prior for Pexpo
inits <- list(
  list(int=0.0, betaquest=0.0, Pexpo=0.3),            
  list(int=1.0, betaquest=1.0, Pexpo=0.5),
  list(int=-1.0, betaquest=-1.0, Pexpo=0.7)
)


# Run the Bayesian analysis -----------------------------------------------
library(R2OpenBUGS)
#Set number of iterations
niterations <- 5000
#Run the Bayesian model
bug.out.exp_non_dif1 <- bugs(data=cyst, 
                             inits=inits, 
                             parameters.to.save=c("int", "betaquest", "ORquest"),
                             n.iter=niterations, 
                             n.burnin=0, 
                             n.thin=1,
                             n.chains=3, 
                             model.file="model_exp_non_dif1.txt", 
                             debug=TRUE,
                             DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.exp_non_dif1)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)


# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_exp_non_dif1 <- data.frame(rbind(bug.out.exp_non_dif1$sims.array[k,1,],bug.out.exp_non_dif1$sims.array[k,2,],bug.out.exp_non_dif1$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

##################################################################
### Compute median estimates with 2.5th and 97.5th percentiles ###
##################################################################
