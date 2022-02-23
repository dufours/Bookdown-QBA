
# Uploading and exploring the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

#Checking number of variables and observations
dim(cyst)

# Checking the first 6 lines
head(cyst)  

# Exploring the quest and test variables
table(cyst$test)

table(cyst$quest)

# Creating a 2 by 2 table
tab <- table(cyst$quest,cyst$test)
# Reordering the values in the 2x2 table to fit epi.2by2() expected data presentation (D+E+, D+E-, D-E+, D-E-)
tab2<- tab[c(2,1),c(2,1)]
# Checking crude associations
library(epiR)
epi.2by2(tab2, method="cohort.count")


# Run a conventional logistic regression model ----------------------------

# Need to create R objects that will be used in the model text file to describe the priors
options(scipen=999) # To avoid use of scientific notation by R

##################################################################################
##################################################################################
# Here you will need to create R objects that will be used to specify the mean and 
# the inverse variance of the prior distributions for your two unknown parameters 
# (i.e. the intercept and the coefficient of the predictor 'quest').
# You also need a R object which will be the number of observations.
##################################################################################
##################################################################################

mu_int <- #TO COMPLETE#                   #mu for the intercept prior distribution    
inv_var_int <- #TO COMPLETE#              #Inverse variance for the intercept prior distribution 
mu_betaquest <- #TO COMPLETE#             #mu for the coefficient prior distribution 
inv_var_betaquest <- #TO COMPLETE#        #Inverse variance for the coefficient prior distribution 

a <- dim(cyst)
num_obs <- a[1] 
num_obs

##################################################################################
##################################################################################
# Here create a text file where the model is written using paste0() and using the 
# R objects previously created.
##################################################################################
##################################################################################

model_vague <- paste0(   #TO COMPLETE#
  
                )

write.table(model_vague, file="model_vague.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

##################################################################################
##################################################################################
# Now you need a number of sets of initial values corresponding to the number of 
# Markov chains that will be ran. Use three chains for now.
# I provided values for a first chain, provide values for the other two chains.
##################################################################################
##################################################################################
inits <- list(
  list(int=0.0, betaquest=0.0),
  list(#TO COMPLETE#
    ),
  list(#TO COMPLETE#
    )
)

# We have everything, we can now run the model using R2OpenBUGS
library(R2OpenBUGS)

#########################
#########################
#Set number of iterations
#########################
#########################

niterations <- #TO COMPLETE#

#Run the Bayesian model
bug.out.vague <- bugs(data=cyst, 
                inits=inits, 
                parameters.to.save=c("int", "betaquest", "ORquest"),
                n.iter=niterations, 
                n.burnin=0, 
                n.thin=1,
                n.chains=3, 
                model.file="model_vague.txt", 
                debug=TRUE,
                DIC=FALSE)

###################################################################################################################
###################################################################################################################
# When you ran this, the OpenBUGS software showed up and ran the analysis. 
# You can close OpenBUGS when its done.
# Check you Environment panel (upper right pannel in RStudio), there is now a R object called bug.out.vague (or any other name you may have used)
# This new R object contains all sort of useful stuff!!!
###################################################################################################################
###################################################################################################################


# MCMC diagnostic ---------------------------------------------------------

#####################################################
#####################################################
# Used the bug.out.vague R object to generate the 
# diagnostic plots.
#####################################################
#####################################################

library(mcmcplots)
#TO COMPLETE#
#TO COMPLETE#


# Applying burn in ---------------------------------

#####################################################
#####################################################
# Set up a burn in period, 
# discard values generated by the MCMC during this period, 
# and combine the post burn in values of the Markov chains.
#####################################################
#####################################################

burnin <- #TO COMPLETE#

# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_conv_vague <- data.frame(rbind(bug.out$sims.array[k,1,],bug.out$sims.array[k,2,],bug.out$sims.array[k,3,]))

############################################################################################################################################
############################################################################################################################################
# Now you have a new R object named estimate_conv_vague (or any other name you may have used) in the Environment panel (upper right pannel in RStudio).
# Check it out, it's the values sampled at each iteration during the MCMC, after the burn in period.
############################################################################################################################################
############################################################################################################################################


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th, and 97.5th)
res_conv_vague <- round(
  t(
    apply(X=estimates_conv_logistic, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_conv_vague

# You could also plot the posterior distributions (here I did it for the predictor coefficient), though it is not mandatory at all. 
# The figure will show up in the Plots panel (lower right panel in RStudio )
std <- sqrt(1/inv_var_betaquest)
plot(density(x=estimates_conv_vague$betaquest), 
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
