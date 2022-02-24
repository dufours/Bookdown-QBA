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
# Rather than proposing prior distributions, I will use exact values (deterministic approach)
Se1 <- 0.70             #VALUE FOR SENSITIVITY IN HEALTHY
Sp1 <- 0.85		        #VALUE FOR SPECIFICITY IN HEALTHY
Se2 <- 0.90             #VALUE FOR SENSITIVITY IN DISEASED
Sp2 <- 0.95		        #VALUE FOR SPECIFICITY IN DISEASED

# Specify the model:
model_exp_dif1 <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							
							#RESPONSE PART
							test[i] ~ dbern(P[i])
							logit(P[i]) <- int + betaquest*expo_true[i] 

              #CREATE VARIABLE INDICATING ACCURACY PARAMETERS TO PICK AS FUNCTION OF OUTCOME TEST
              diff[i] <- test[i]+1
							
							#MEASUREMENT ERROR PART
							expo_true[i] ~ dbern(Pexpo)
							quest[i] ~ dbern(P_obs_expo[i])
							P_obs_expo[i] <- Se[diff[i]]*expo_true[i] + (1-Sp[diff[i]])*(1-expo_true[i])
					
			}

#===   EXTRA CALCULATION TO GET OR ===#
			ORquest <- exp(betaquest)
			
#===	PRIOR	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		#PRIOR INTERCEPT
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#PRIOR COEFFICIENT
			Pexpo ~ dbeta(1.0, 1.0)                                     #FLAT PRIOR ON PREVALENCE OF TRUE EXPOSURE
			
			Se[1] <- ", Se1, "		    #PRIOR FOR SENSITIVITY HEALTHY
			Sp[1] <- ", Sp1, "		    #PRIOR FOR SPECIFICITY HEALTHY
			Se[2] <- ", Se2, "		    #PRIOR FOR SENSITIVITY DISEASED
			Sp[2] <- ", Sp2, "		    #PRIOR FOR SPECIFICITY DISEASED
			
			}")

#write to temporary text file
write.table(model_exp_dif1, file="model_exp_dif1.txt", quote=FALSE, sep="", row.names=FALSE,
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
bug.out.exp_dif1 <- bugs(data=cyst, 
                             inits=inits, 
                             parameters.to.save=c("int", "betaquest", "ORquest"),
                             n.iter=niterations, 
                             n.burnin=0, 
                             n.thin=1,
                             n.chains=3, 
                             model.file="model_exp_dif1.txt", 
                             debug=TRUE,
                             DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.exp_dif1)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_exp_dif1 <- data.frame(rbind(bug.out.exp_dif1$sims.array[k,1,],bug.out.exp_dif1$sims.array[k,2,],bug.out.exp_dif1$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_exp_dif1 <- round(
  t(
    apply(X=estimates_exp_dif1, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_exp_dif1
