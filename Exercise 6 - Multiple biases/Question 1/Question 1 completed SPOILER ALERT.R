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
#Sensitivity and specificity of the test used to diagnose exposure
Se <- 0.80             #VALUE FOR SENSITIVITY 
Sp <- 0.90		        #VALUE FOR SPECIFICITY 

# Selection probabilities in exposed and unexposed individuals
S_dis_exp <- 0.75             #VALUE FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
S_dis_unexp <- 0.25		        #VALUE FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
S_undis_exp <- 0.30		        #VALUE FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
S_undis_unexp <- 0.10		      #VALUE FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS 

# Association between confounder and disease + prevalence of confounder in exposed and unexposed

OR_ud <- 5         #Association between confounder and disease
PU_exp <- 0.60     # Prevalence of confounder in exposed
PU_unexp <- 0.20   # Prevalence of confounder in unexposed
  

# Specify the model:
model_mult <- paste0("model{

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

#===  EXTRA CALCULATION TO GET THE ADJUSTED OR ===#
		
			#OR ADJUSTED FOR EXPOSURE MISCLASSIFICATION
			OR_adj_exp_mis <- exp(betaquest)
			
			#OR ADJUSTED FOR SELECTION BIAS AND EXPOSURE MISCLASSIFICATION
			OR_sel <- (S_dis_unexp*S_undis_exp)/(S_dis_exp*S_undis_unexp)
			OR_adj_exp_sel <- OR_adj_exp_mis*OR_sel
			
			#OR ADJUSTED FOR UNMEASURED CONFOUNDER, SELECTION BIAS AND EXPOSURE MISCLASSIFICATION
			OR_adj_conf_exp_sel <- OR_adj_exp_sel*( (OR_ud*PU_unexp+1-PU_unexp)/(OR_ud*PU_exp+1-PU_exp) )
			
#===	PRIOR	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		#PRIOR INTERCEPT
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#PRIOR COEFFICIENT
			Pexpo ~ dbeta(1.0, 1.0)                                     #FLAT PRIOR ON PREVALENCE OF TRUE EXPOSURE
			
			Se <- ", Se, "		    #VALUE FOR SENSITIVITY
			Sp <- ", Sp, "		    #VALUE FOR SPECIFICITY
			
			S_dis_exp <- ", S_dis_exp, "             #VALUE FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
      S_dis_unexp <- ", S_dis_unexp, "	        #VALUE FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
      S_undis_exp <- ", S_undis_exp, "	        #VALUE FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
      S_undis_unexp <- ", S_undis_unexp, "		      #VALUE FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS
      
      OR_ud <- ", OR_ud, "         #Association between confounder and disease
      PU_exp <- ", PU_exp, "     # Prevalence of confounder in exposed
      PU_unexp <- ", PU_unexp, "   # Prevalence of confounder in unexposed
			
			}")

#write to temporary text file
write.table(model_mult, file="model_mult.txt", quote=FALSE, sep="", row.names=FALSE,
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
bug.out.mult <- bugs(data=cyst, 
                             inits=inits, 
                             parameters.to.save=c("int", "betaquest", "OR_adj_exp_mis", "OR_adj_exp_sel", "OR_adj_conf_exp_sel"),
                             n.iter=niterations, 
                             n.burnin=0, 
                             n.thin=1,
                             n.chains=3, 
                             model.file="model_mult.txt", 
                             debug=TRUE,
                             DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.mult)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_mult <- data.frame(rbind(bug.out.mult$sims.array[k,1,],bug.out.mult$sims.array[k,2,],bug.out.mult$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_mult <- round(
  t(
    apply(X=estimates_mult, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_mult

# Plot fully adjusted OR 
plot(density(x=estimates_mult$OR_adj_conf_exp_sel), 
     main="Adjusted (deterministic) vs. biased OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 15),
     xlim=c(0, 1.0)
)
# Posterior distribution of exposure and selection adjusted OR  
lines(density(estimates_mult$OR_adj_exp_sel), col="blue")

# Posterior distribution of exposure adjusted OR  
lines(density(estimates_mult$OR_adj_exp_mis), col="purple")

# Posterior distribution of biased OR  (NEED TO RUN Q1 in exercise 1 BEFOREHAND)
lines(density(estimates_conv_vague$ORquest), col="red", lty="dashed")

