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


# Priors for the bias parameters ------------------------------------------
# Now I need to find corresponding prior distributions  
library(epiR) 
#Se
rval.se <- epi.betabuster(mode=0.80, conf=0.95, greaterthan=T, x=0.75)  
rval.se$shape1                 
rval.se$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.se$shape1, shape2=rval.se$shape2), from=0, to=1, 
      main="Prior for test sensitivity", xlab = "Sensitivity", ylab = "Density")

#Sp
rval.sp <- epi.betabuster(mode=0.90, conf=0.95, greaterthan=T, x=0.85)  
rval.sp$shape1                 
rval.sp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.sp$shape1, shape2=rval.sp$shape2), from=0, to=1, 
      main="Prior for test specificity", xlab = "Specificity", ylab = "Density")

# Specify the model:
model_exp_non_dif2 <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							
							#RESPONSE PART
							test[i] ~ dbern(P[i])
							logit(P[i]) <- int + betaquest*expo_true[i] 
							
							#MEASUREMENT ERROR PART
							expo_true[i] ~ dbern(Pexpo)
							quest[i] ~ dbern(P_obs_expo[i])
							P_obs_expo[i] <- Se*expo_true[i] + (1-Sp)*(1-expo_true[i])
					
			}

#===   EXTRA CALCULATION TO GET OR ===#
			ORquest <- exp(betaquest)
			
#===	PRIOR	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		#PRIOR INTERCEPT
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#PRIOR COEFFICIENT
			Pexpo ~ dbeta(1.0, 1.0)                                     #FLAT PRIOR ON PREVALENCE OF TRUE EXPOSURE
			
			Se ~ dbeta(", rval.se$shape1,",",rval.se$shape2, ")		    #PRIOR FOR SENSITIVITY
			Sp ~ dbeta(", rval.sp$shape1,",",rval.sp$shape2, ")		    #PRIOR FOR SPECIFICITY
                             
      }")

#write to temporary text file
write.table(model_exp_non_dif2, file="model_exp_non_dif2.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

# Generate initial values -------------------------------------------------

# In this case, I do not need initial values for the bias parameters, since we specified exact values for them (vs distributions)
# But I do need a prior for Pexpo
inits <- list(
  list(int=0.0, betaquest=0.0, Pexpo=0.3, Se=0.8, Sp=0.9),            
  list(int=1.0, betaquest=1.0, Pexpo=0.5, Se=0.7, Sp=0.8),
  list(int=-1.0, betaquest=-1.0, Pexpo=0.7, Se=0.9, Sp=0.8)
)


# Run the Bayesian analysis -----------------------------------------------
library(R2OpenBUGS)
#Set number of iterations
niterations <- 15000
#Run the Bayesian model
bug.out.exp_non_dif2 <- bugs(data=cyst, 
                             inits=inits, 
                             parameters.to.save=c("int", "betaquest", "ORquest"),
                             n.iter=niterations, 
                             n.burnin=0, 
                             n.thin=1,
                             n.chains=3, 
                             model.file="model_exp_non_dif2.txt", 
                             debug=TRUE,
                             DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.exp_non_dif2)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_exp_non_dif2 <- data.frame(rbind(bug.out.exp_non_dif2$sims.array[k,1,],bug.out.exp_non_dif2$sims.array[k,2,],bug.out.exp_non_dif2$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_exp_non_dif2 <- round(
  t(
    apply(X=estimates_exp_non_dif2, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_exp_non_dif2

# Plot unadjusted and adjusted OR (NEED TO RUN Q1 AND Q2 BEFOREHAND)
plot(density(x=estimates_exp_non_dif2$ORquest), 
     main="Adj. stochastic (black) vs. adj. deterministic (blue) OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 7),
     xlim=c(0, 2.0)
)
# Posterior distribution of biased OR - deterministic    #####NEED TO RUN QUESTION 2 FIRST####
lines(density(estimates_exp_non_dif1$ORquest), col="blue")

