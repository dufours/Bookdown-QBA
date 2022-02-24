# Uploading and exploring the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

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
# Rather than proposing prior distributions, I will use exact values (deterministic approach)
S_dis_exp <- 0.75             #VALUE FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
S_dis_unexp <- 0.75		        #VALUE FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
S_undis_exp <- 0.30		        #VALUE FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
S_undis_unexp <- 0.10		      #VALUE FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS 


# Specify the model ------------------------------------------------------

model_sel1 <- paste0("model{

#===	LIKELIHOOD	===#
			for( i in 1 : ",num_obs," ) {
							test[i] ~ dbern(p[i])
							logit(p[i]) <- int + betaquest*quest[i]  
			}

#===  EXTRA CALCULATION TO COMPUTE OR_sel ===#
			OR_sel <- (S_dis_unexp*S_undis_exp)/(S_dis_exp*S_undis_unexp)
			
#===  EXTRA CALCULATION TO GET OR_biased and OR_adj ===#
			OR_biased <- exp(betaquest)
			OR_adj <- OR_biased*OR_sel
			
#===	PRIORS	===#		
			int ~ dnorm(",mu_int,", ",inv_var_int,")		                #FLAT PRIOR
			betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")	#FLAT PRIOR
		
			S_dis_exp <- ", S_dis_exp,"		              #VALUE FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
			S_dis_unexp <- ", S_dis_unexp,"		          #VALUE FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
			S_undis_exp <- ", S_undis_exp,"		            #VALUE FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
			S_undis_unexp <- ", S_undis_unexp,"		          #VALUE FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS 

		}")

#write to temporary text file
write.table(model_sel1, file="model_sel1.txt", quote=FALSE, sep="", row.names=FALSE,
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
#Run the Bayesian model
bug.out.sel2 <- bugs(data=cyst, 
                    inits=inits, 
                    parameters.to.save=c("int", "betaquest", "OR_biased", "OR_adj"),
                    n.iter=niterations, 
                    n.burnin=0, 
                    n.thin=1,
                    n.chains=3, 
                    model.file="model_sel1.txt", 
                    debug=TRUE,
                    DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.sel2)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_sel2 <- data.frame(rbind(bug.out.sel2$sims.array[k,1,],bug.out.sel2$sims.array[k,2,],bug.out.sel2$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_sel2 <- round(
  t(
    apply(X=estimates_sel2, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_sel2

# Plot the posterior distributions for the adjusted vs the biased OR (red)
plot(density(x=estimates_sel2$OR_adj), 
     main="Adjusted (black) vs. biased (red dashed) OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 7),
     xlim=c(0, 2.5)
)
# Posterior distribution of OR_biased
lines(density(estimates_sel2$OR_biased), col="red", lty="dashed")



#We could plot the log(odds) to better appraise the differences
#Need to compute the adjusted OR on log(odds) scale first
estimates_sel2$beta_adj <- log(estimates_sel2$OR_adj)
# Plot the posterior distributions for the adjusted vs the biased OR (red)
plot(density(x=estimates_sel2$beta_adj), 
     main="Adjusted (black) vs. biased (red dashed) log(OR)",
     xlab="OR", ylab="Density",
     ylim=c(0, 4),
     xlim=c(-1.5, 1)
)
# Posterior distribution of OR_biased
lines(density(estimates_sel2$betaquest), col="red", lty="dashed")
