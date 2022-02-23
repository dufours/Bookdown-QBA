# Uploading and exploring the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

# Creating a 2 by 2 table
tab <- table(cyst$quest,cyst$test)
# Reordering the values in the 2x2 table to fit epi.2by2() expected data presentation (D+E+, D+E-, D-E+, D-E-)
tab2<- tab[c(2,1),c(2,1)]
# Checking crude associations
library(epiR)
epi.2by2(tab2, method="cohort.count")

#Create a "a" object with number of rows and columns
a <- dim(cyst)
#Save and check the first value of the "a" object
num_obs <- a[1]

# Priors for intercept and coefficient -------------------------------------------------------------
options(scipen=999) # To avoid use of scientific notation by R

mu_int <- 0.0          #mu for the intercept    
inv_var_int <- 0.0001  #Inverse variance for the intercept
mu_betaquest <- 0.0    #mu for the coefficient
inv_var_betaquest <- 0.0001 #Inverse variance for the coefficient


# Values for the bias parameters ------------------------------------------
##########################################################################################################
# Rather than proposing prior distributions, you could first use exact values (deterministic approach) ###
##########################################################################################################

###TO COMPLETE###


# Specify the model ------------------------------------------------------

#############################################################
### Specify the model's likelihood function              ###
### Add the values for Se and Sp                          ###
#############################################################

model_mis_non_dif <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							
							###TO COMPLETE###
			}

#===   EXTRA CALCULATION TO GET OR ===#
			ORquest <- exp(betaquest)
			
#===	PRIOR	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		#FLAT PRIOR
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#FLAT PRIOR
			
			###TO COMPLETE###		    #VALUE FOR SENSITIVITY
			###TO COMPLETE###		    #VALUE FOR SPECIFICITY
			}")

#write to temporary text file
write.table(model_mis_non_dif, file="model_mis_non_dif.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

# Generate initial values -------------------------------------------------

# In this case, I do not need initial values for the bias parameters, since we specified exact values for them (vs distributions)
inits <- list(
  list(int=0.0, betaquest=0.0),            
  list(int=1.0, betaquest=1.0),
  list(int=-1.0, betaquest=-1.0)
)


# Run the Bayesian analysis -----------------------------------------------

library(R2OpenBUGS)
#Set number of iterations
niterations <- 5000

####################################################
### Below, specify the parameters to monitor    ####
####################################################

#Run the Bayesian model
bug.out.mis_non_dif1 <- bugs(data=cyst, 
                     inits=inits, 
                     parameters.to.save=c(###TO COMPLETE###),
                     n.iter=niterations, 
                     n.burnin=0, 
                     n.thin=1,
                     n.chains=3, 
                     model.file="model_mis_non_dif.txt", 
                     debug=TRUE,
                     DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.mis_non_dif1)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000

####################################################
### Adjust burn in and thinning as you deem fit  ###
####################################################

burnin <- ###TO COMPLETE###
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=###TO COMPLETE###)
# Combine chains post-burn in
estimates_mis_non_dif1 <- data.frame(rbind(bug.out.mis_non_dif1$sims.array[k,1,],bug.out.mis_non_dif1$sims.array[k,2,],bug.out.mis_non_dif1$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_mis_non_dif1 <- round(
  t(
    apply(X=estimates_mis_non_dif1, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_mis_non_dif1


