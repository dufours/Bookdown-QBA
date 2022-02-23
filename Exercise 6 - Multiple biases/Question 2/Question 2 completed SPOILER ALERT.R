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

#Sensitivity and specificity of the test used to diagnose exposure
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

# Selection probabilities in exposed and unexposed individuals
#S(Dis+ Exp+)
rval.dis.exp <- epi.betabuster(mode=0.75, conf=0.95, greaterthan=T, x=0.70)  
rval.dis.exp$shape1                 
rval.dis.exp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.dis.exp$shape1, shape2=rval.dis.exp$shape2), from=0, to=1, 
      main="Prior for S(d+e+)", xlab = "Proportion", ylab = "Density")

#S(Dis+ Exp-)
rval.dis.unexp <- epi.betabuster(mode=0.25, conf=0.95, greaterthan=F, x=0.30)  
rval.dis.unexp$shape1                 
rval.dis.unexp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.dis.unexp$shape1, shape2=rval.dis.unexp$shape2), from=0, to=1, 
      main="Prior for S(d+e-)", xlab = "Proportion", ylab = "Density")

#S(Dis- Exp+)
rval.undis.exp <- epi.betabuster(mode=0.30, conf=0.95, greaterthan=F, x=0.35)  
rval.undis.exp$shape1                 
rval.undis.exp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.undis.exp$shape1, shape2=rval.undis.exp$shape2), from=0, to=1, 
      main="Prior for S(d-e+)", xlab = "Proportion", ylab = "Density")

#S(Dis- Exp-)
rval.undis.unexp <- epi.betabuster(mode=0.10, conf=0.95, greaterthan=F, x=0.15)  
rval.undis.unexp$shape1                 
rval.undis.unexp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.undis.unexp$shape1, shape2=rval.undis.unexp$shape2), from=0, to=1, 
      main="Prior for S(d-e-)", xlab = "Proportion", ylab = "Density")

# Association between confounder and disease + prevalence of confounder in exposed and unexposed

#1) association
mu_log_OR_ud <- log(5)         #Mean association (on log odds scale) between confounder and disease (corresponding to an OR of 5.0)   
inv_var_log_OR_ud <- 1/0.2     #Inverse variance for the association (on log odds scale) between confounder and disease (corresponding to a variance of 0.2)

#2) I need priors for the prevalence of the confounder in the exposed
library(epiR) 
rval1 <- epi.betabuster(mode=0.60, conf=0.95, greaterthan=T, x=0.50)  
rval1$shape1                #View the a shape parameter 
rval1$shape2                ##View the b shape parameter
#plot the prior distribution
curve(dbeta(x, shape1=rval1$shape1, shape2=rval1$shape2), from=0, to=1, 
      main="Prior for prevalence of confounder in exposed", xlab = "Proportion", ylab = "Density")

#3) I need priors for the prevalence of the confounder in the unexposed
rval2 <- epi.betabuster(mode=0.20, conf=0.95, greaterthan=F, x=0.30)  
rval2$shape1                #View the a shape parameter 
rval2$shape2                ##View the b shape parameter 
#plot the prior distribution
curve(dbeta(x, shape1=rval2$shape1, shape2=rval2$shape2), from=0, to=1, 
      main="Prior for prevalence of confounder in unexposed", xlab = "Proportion", ylab = "Density")

# Specify the model:
model_mult2 <- paste0("model{

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
			
			Se ~ dbeta(", rval.se$shape1,",",rval.se$shape2, ")		    #PRIOR FOR SENSITIVITY
			Sp ~ dbeta(", rval.sp$shape1,",",rval.sp$shape2, ")		    #PRIOR FOR SPECIFICITY

			S_dis_exp ~ dbeta(", rval.dis.exp$shape1,",",rval.dis.exp$shape2, ")		       #PRIOR FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
			S_dis_unexp ~ dbeta(", rval.dis.unexp$shape1,",",  rval.dis.unexp$shape2, ")		  #PRIOR FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
			S_undis_exp ~ dbeta(",  rval.undis.exp$shape1,",", rval.undis.exp$shape2, ")		     #PRIOR FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
			S_undis_unexp ~ dbeta(", rval.undis.unexp$shape1,",", rval.undis.unexp$shape2, ")		                #PRIOR FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS 
      
			log_OR_ud ~ dnorm(", mu_log_OR_ud,"," , inv_var_log_OR_ud,")   #PRIOR FOR log_OR_ud		
			OR_ud <- exp(log_OR_ud)	                                       #LINK OR_ud TO THE LOG_OR_UD
			PU_exp ~ dbeta(", rval1$shape1,",", rval1$shape2, ")		       #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE EXPOSED INDIVIDUALS 
			PU_unexp ~ dbeta(", rval2$shape1,",", rval2$shape2, ")		     #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE UNEXPOSED INDIVIDUALS 	
			
			}")

#write to temporary text file
write.table(model_mult2, file="model_mult2.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

# Generate initial values -------------------------------------------------

# In this case, I do need initial values for the bias parameters
inits <- list(
  list(int=0.0, betaquest=0.0, Pexpo=0.3, Se=0.8, Sp=0.9, S_dis_exp=0.25, S_dis_unexp=0.20, S_undis_exp=0.20, S_undis_unexp=0.20, log_OR_ud=0, PU_exp=0.20, PU_unexp=0.20),            
  list(int=1.0, betaquest=1.0, Pexpo=0.5, Se=0.7, Sp=0.8, S_dis_exp=0.35, S_dis_unexp=0.30, S_undis_exp=0.30, S_undis_unexp=0.30, log_OR_ud=1, PU_exp=0.30, PU_unexp=0.30),
  list(int=-1.0, betaquest=-1.0, Pexpo=0.7, Se=0.9, Sp=0.8, S_dis_exp=0.45, S_dis_unexp=0.40, S_undis_exp=0.40, S_undis_unexp=0.40, log_OR_ud=-1, PU_exp=0.10, PU_unexp=0.10)
)


# Run the Bayesian analysis -----------------------------------------------
library(R2OpenBUGS)
#Set number of iterations
niterations <- 15000
#Run the Bayesian model
bug.out.mult2 <- bugs(data=cyst, 
                             inits=inits, 
                             parameters.to.save=c("int", "betaquest", "OR_adj_exp_mis", "OR_adj_exp_sel", "OR_adj_conf_exp_sel"),
                             n.iter=niterations, 
                             n.burnin=0, 
                             n.thin=1,
                             n.chains=3, 
                             model.file="model_mult2.txt", 
                             debug=TRUE,
                             DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.mult2)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_mult2 <- data.frame(rbind(bug.out.mult2$sims.array[k,1,],bug.out.mult2$sims.array[k,2,],bug.out.mult2$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_mult2 <- round(
  t(
    apply(X=estimates_mult2, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_mult2

# Plot fully adjusted OR 
plot(density(x=estimates_mult2$OR_adj_conf_exp_sel), 
     main="Adjusted (deterministic) vs. biased OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 8),
     xlim=c(0, 1.0)
)
# Posterior distribution of exposure and selection adjusted OR  
lines(density(estimates_mult2$OR_adj_exp_sel), col="blue")

# Posterior distribution of exposure adjusted OR  
lines(density(estimates_mult2$OR_adj_exp_mis), col="purple")

# Posterior distribution of biased OR  (NEED TO RUN Q1 in exercise 1 BEFOREHAND)
lines(density(estimates_conv_vague$ORquest), col="red", lty="dashed")

