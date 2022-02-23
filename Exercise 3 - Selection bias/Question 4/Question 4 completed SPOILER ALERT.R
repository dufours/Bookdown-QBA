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


# Priors for the bias parameters ------------------------------------------
library(epiR) 
#S(Dis+ Exp+)
rval.dis.exp <- epi.betabuster(mode=0.75, conf=0.95, greaterthan=T, x=0.55)  
rval.dis.exp$shape1                 
rval.dis.exp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.dis.exp$shape1, shape2=rval.dis.exp$shape2), from=0, to=1, 
      main="Prior for S(d+e+)", xlab = "Proportion", ylab = "Density")

#S(Dis+ Exp-)
rval.dis.unexp <- epi.betabuster(mode=0.75, conf=0.95, greaterthan=T, x=0.55)  
rval.dis.unexp$shape1                 
rval.dis.unexp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.dis.unexp$shape1, shape2=rval.dis.unexp$shape2), from=0, to=1, 
      main="Prior for S(d+e-)", xlab = "Proportion", ylab = "Density")

#S(Dis- Exp+)
rval.undis.exp <- epi.betabuster(mode=0.30, conf=0.95, greaterthan=F, x=0.50)  
rval.undis.exp$shape1                 
rval.undis.exp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.undis.exp$shape1, shape2=rval.undis.exp$shape2), from=0, to=1, 
      main="Prior for S(d-e+)", xlab = "Proportion", ylab = "Density")

#S(Dis- Exp-)
rval.undis.unexp <- epi.betabuster(mode=0.10, conf=0.95, greaterthan=F, x=0.30)  
rval.undis.unexp$shape1                 
rval.undis.unexp$shape2                   
#plot the prior distribution
curve(dbeta(x, shape1=rval.undis.unexp$shape1, shape2=rval.undis.unexp$shape2), from=0, to=1, 
      main="Prior for S(d-e-)", xlab = "Proportion", ylab = "Density")


# Specify the model -------------------------------------------------------

# Specify the model:
model_sel4 <- paste0("model{

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
		
			S_dis_exp ~ dbeta(", rval.dis.exp$shape1,",",rval.dis.exp$shape2, ")		       #PRIOR FOR SELECTION PROBABILITY IN D+E+ INDIVIDUALS 
			S_dis_unexp ~ dbeta(", rval.dis.unexp$shape1,",",  rval.dis.unexp$shape2, ")		  #PRIOR FOR SELECTION PROBABILITY IN D+E- INDIVIDUALS  
			S_undis_exp ~ dbeta(",  rval.undis.exp$shape1,",", rval.undis.exp$shape2, ")		     #PRIOR FOR SELECTION PROBABILITY IN D-E+ INDIVIDUALS 
			S_undis_unexp ~ dbeta(", rval.undis.unexp$shape1,",", rval.undis.unexp$shape2, ")		                #PRIOR FOR SELECTION PROBABILITY IN D-E- INDIVIDUALS 

		}")

#write to temporary text file
write.table(model_sel4, file="model_sel4.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)


# Specify initial values --------------------------------------------------

inits <- list(
  list(int=0.0, betaquest=0.0, S_dis_exp=0.20, S_dis_unexp=0.20, S_undis_exp=0.20, S_undis_unexp=0.20),            
  list(int=1.0, betaquest=1.0, S_dis_exp=0.30, S_dis_unexp=0.30, S_undis_exp=0.30, S_undis_unexp=0.30),
  list(int=-1.0, betaquest=-1.0, S_dis_exp=0.40, S_dis_unexp=0.40, S_undis_exp=0.40, S_undis_unexp=0.40)
)


# Run the model -----------------------------------------------------------

library(R2OpenBUGS)
#Set number of iterations
niterations <- 5000
#Run the Bayesian model
bug.out.sel4 <- bugs(data=cyst, 
                    inits=inits, 
                    parameters.to.save=c("int", "betaquest", "OR_sel", "OR_biased", "OR_adj"),
                    n.iter=niterations, 
                    n.burnin=0, 
                    n.thin=1,
                    n.chains=3, 
                    model.file="model_sel4.txt", 
                    debug=TRUE,
                    DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.sel4)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1)
# Combine chains post-burn in
estimates_sel4 <- data.frame(rbind(bug.out.sel4$sims.array[k,1,],bug.out.sel4$sims.array[k,2,],bug.out.sel4$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_sel4 <- round(
  t(
    apply(X=estimates_sel4, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_sel4

# Plot the posterior distributions for the adjusted vs the biased OR (red)
plot(density(x=estimates_sel4$OR_adj), 
     main="Adj. stochastic wide (black) vs. narrow (blue) priors vs. adj. deterministic (green) vs. biased (red dashed) OR",
     xlab="OR", ylab="Density",
     ylim=c(0, 7),
     xlim=c(0, 2.0)
)
# Posterior distribution of OR_adj    #####NEED TO RUN QUESTION 3 FIRST####
lines(density(estimates_sel3$OR_adj), col="blue")
# Posterior distribution of biased OR - deterministic    #####NEED TO RUN QUESTION 2 FIRST####
lines(density(estimates_sel2$OR_adj), col="green")
# Posterior distribution of biased OR
lines(density(estimates_sel4$OR_biased), col="red", lty="dashed")

