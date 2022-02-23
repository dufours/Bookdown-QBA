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


##############################################################
### Below, I tried to use informative prior distribution   ###
##############################################################

#mu_int <- 0.0          #mu for the intercept    
#inv_var_int <- 0.0001  #Inverse variance for the intercept
#mu_betaquest <- 0.0    #mu for the coefficient
#inv_var_betaquest <- 0.0001 #Inverse variance for the coefficient

int.low <- -2.0
int.hi <- 2.0
#PRIOR = PREVALENCE IN UNEXPOSED >0.12 and <0.88
beta.low <- -2.0
beta.hi <- 2.0
#PRIOR = OR BETWEEN 0.14 AND 7.4

# Values for the bias parameters ------------------------------------------
# Now I need to find corresponding prior distributions  
library(epiR) 
#Se
rval.se <- epi.betabuster(mode=0.90, conf=0.95, greaterthan=T, x=0.80)  
rval.se$shape1                 
rval.se$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.se$shape1, shape2=rval.se$shape2), from=0, to=1, 
main="Prior for test sensitivity", xlab = "Sensitivity", ylab = "Density")

#Sp
rval.sp <- epi.betabuster(mode=0.97, conf=0.95, greaterthan=T, x=0.90)  
rval.sp$shape1                 
rval.sp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.sp$shape1, shape2=rval.sp$shape2), from=0, to=1, 
main="Prior for test specificity", xlab = "Specificity", ylab = "Density")


# Specify the model ------------------------------------------------------

model_mis_non_dif <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							test[i] ~ dbern(P_obs[i])
							P_obs[i] <- P_true[i]*Se+(1-P_true[i])*(1-Sp)	#MEASUREMENT ERROR PART
							logit(P_true[i]) <- int + betaquest*quest[i]  #RESPONSE PART
			}

#===   EXTRA CALCULATION TO GET OR ===#
			ORquest <- exp(betaquest)
			
#===	PRIOR	===#		
			int ~ dunif(",int.low,", ",int.hi,")		#PRIOR INTERCEPT 
			betaquest ~ dunif(",beta.low,", ",beta.hi,")	#PRIOR COEFFICIENT
			
			Se ~ dbeta(", rval.se$shape1,",",rval.se$shape2, ")		    #PRIOR FOR SENSITIVITY
			Sp ~ dbeta(", rval.sp$shape1,",",rval.sp$shape2, ")		    #PRIOR FOR SPECIFICITY
			}")

#write to temporary text file
write.table(model_mis_non_dif, file="model_mis_non_dif.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

# Generate initial values -------------------------------------------------

# In this case, I do need initial values for the bias parameters also
inits <- list(
  list(int=0.0, betaquest=0.0, Se=0.8, Sp=0.8),            
  list(int=1.0, betaquest=1.0, Se=0.9, Sp=0.9),
  list(int=-1.0, betaquest=-1.0, Se=0.7, Sp=0.7)
)


# Run the Bayesian analysis -----------------------------------------------

library(R2OpenBUGS)

####################################################################
### Below I had to finetune a lot the estimation process         ###
### I used 15,000 iterations, and a burn in of 5000              ###
####################################################################

#Set number of iterations
niterations <- 15000
#Run the Bayesian model
bug.out.mis_non_dif3 <- bugs(data=cyst, 
                     inits=inits, 
                     parameters.to.save=c("int", "betaquest", "ORquest"),
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
bug.mcmc <- as.mcmc(bug.out.mis_non_dif3)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------
# Set burn-in to 5000
burnin <- 5000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_mis_non_dif3 <- data.frame(rbind(bug.out.mis_non_dif3$sims.array[k,1,],bug.out.mis_non_dif3$sims.array[k,2,],bug.out.mis_non_dif3$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_mis_non_dif3 <- round(
  t(
    apply(X=estimates_mis_non_dif3, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_mis_non_dif3

# Plot unadjusted and adjusted OR (NEED TO RUN Q1 AND Q2 BEFOREHAND)
plot(density(x=estimates_mis_non_dif3$ORquest), 
     main="Adj. stochastic (black) vs. adj. deterministic (blue) vs. biased (red dashed) OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 7),
     xlim=c(0, 2.0)
)
# Posterior distribution of biased OR - deterministic    #####NEED TO RUN QUESTION 2 FIRST####
lines(density(estimates_mis_non_dif2$ORquest), col="blue")
# Posterior distribution of biased OR
lines(density(estimates_mis_non_dif1$ORquest), col="red", lty="dashed")


