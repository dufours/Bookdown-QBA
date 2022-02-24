# Uploading the dataset ---------------------------------------------------

cyst <- read.csv(file="Data/cyst_bias_study.csv", header=TRUE, sep=";")

#I need a R object that will be used in the model text file for the number of observations
a <- dim(cyst)
num_obs <- a[1] 

# Priors for intercept an coefficient of quest ----------------------------

#I need to create R objects that will be used in the model text file to describe the priors
options(scipen=999)           #To avoid use of scientific notation by R
mu_int <- 0.0                 #mu for the intercept    
inv_var_int <- 0.0001         #Inverse variance for the intercept
mu_betaquest <- 0.0           #mu for the coefficient
inv_var_betaquest <- 0.0001   #Inverse variance for the coefficient


# Priors for bias parameters ----------------------------------------------

#1) I need priors for the log of the unmeasured confounder-disease OR

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

#3) I need priors for the prevalence of the confounder in the exposed
rval2 <- epi.betabuster(mode=0.20, conf=0.95, greaterthan=F, x=0.30)  
rval2$shape1                #View the a shape parameter 
rval2$shape2                ##View the b shape parameter 
#plot the prior distribution
curve(dbeta(x, shape1=rval2$shape1, shape2=rval2$shape2), from=0, to=1, 
      main="Prior for prevalence of confounder in unexposed", xlab = "Proportion", ylab = "Density")


# Run Bayesian model using Lash's adjustment for unmeasured confounder ----------------------------

#Generate the .txt file with the model
# Specify the model:
model_unm_conf_McCandless <- paste0("model{

#===    LIKELIHOOD  ===#
            for( i in 1 : ",num_obs," ) {
                            test[i] ~ dbern(p[i])
                            logit(p[i]) <- int + betaquest*quest[i] + betaUC*UC[i]      #UC IS AN UNMEASURED CONFOUNDER (LATENT)
                            
                            UC[i] ~ dbern(Pu[i])           #THE VALUE OF UC HAS PROBABILITY Pu OF BEING 1
                            logit(Pu[i]) <- gamma + gammaquest*quest[i]      #Pu VARIES AS FUNCTION OF THE EXPOSURE STATUS
            }
            
            gamma <- log(1/PU_unexp - 1)*-1             #LINK BETWEEN gamma and PU_unexp
            gammaquest <-  (log(1/PU_exp - 1)*-1)-gamma        #LINK BETWEEN gammaquest AND PU_exp


#===    PRIORS  ===#        
            int ~ dnorm(",mu_int,", ",inv_var_int,")                        #FLAT PRIOR
            betaquest ~ dnorm(",mu_betaquest,", ",inv_var_betaquest,")  #FLAT PRIOR
            
            betaUC ~ dnorm(", mu_log_OR_ud,"," , inv_var_log_OR_ud,")   #PRIOR FOR betaUC WHICH IS SIMPLY OR_ud, BUT ON THE LOG SCALE       
            
            PU_exp ~ dbeta(", rval1$shape1,",", rval1$shape2, ")               #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE EXPOSED INDIVIDUALS 
            PU_unexp ~ dbeta(", rval2$shape1,",", rval2$shape2, ")           #PRIOR FOR PREVALENCE OF THE CONFOUNDER IN THE UNEXPOSED INDIVIDUALS   


#===   EXTRA CALCULATION TO GET OR(adj) ===#
            ORquest<-exp(betaquest)
        }")

#write to temporary text file
write.table(model_unm_conf_McCandless, file="model_unm_conf_McCandless.txt", quote=FALSE, sep="", row.names=FALSE,
            col.names=FALSE)

#Now I need a number of sets of initial values corresponding to the number of Markov chains that will be ran. I plan to use 3 chains, thus 3 sets of initial values
#Initializing values for 3 chains
inits <- list(
  list(int=0.0, betaquest=0.0, log_OR_ud=0, PU_exp=0.20, PU_unexp=0.20),            
  list(int=1.0, betaquest=1.0, log_OR_ud=1.0, PU_exp=0.30, PU_unexp=0.30),
  list(int=-1.0, betaquest=-1.0, log_OR_ud=-1.0, PU_exp=0.10, PU_unexp=0.10)
)

# We have everything, we can now run the model using R2OpenBUGS
library(R2OpenBUGS)

########################################################################################################################################
### I first ran the model below with 5000 iterations. Autocorrelation was high, so I  have therefore increased                       ###
### the number of iterations to 15,000, so we still have enough "usable"  iterations per chain.                                      ###                                                            ### 
### IT WILL TAKE A LONG TIME TO RUN!!!                                                                                               ###
########################################################################################################################################

#Set number of iterations
niterations <- 15000
#Run the Bayesian model
bug.out.unm_conf_McCandless <- bugs(data=cyst, 
                      inits=inits, 
                      parameters.to.save=c("int", "betaquest", "ORquest"),  
                      n.iter=niterations, 
                      n.burnin=0, 
                      n.thin=1,                                            
                      n.chains=3, 
                      model.file="model_unm_conf_McCandless.txt", 
                      debug=TRUE,
                      DIC=FALSE)


# MCMC diagnostic ---------------------------------------------------------

# I need to generate the diagnostic plots for the MCMC 
library(mcmcplots)
bug.mcmc <- as.mcmc(bug.out.unm_conf_McCandless)
mcmcplot(bug.mcmc, title="Diagnostic plots")
effectiveSize(bug.mcmc)

# Applying burn in ---------------------------------

# Now we need to set up a burn in period and discard values generated by the MCMC during this period
# Set burn-in to 1000
burnin <- 1000
# Create a sequence of values from 1st post-burn in value to last one.
k <- seq(from=(burnin+1), to=niterations, by=1) 
# Combine chains post-burn in
estimates_unm_conf_McCandless <- data.frame(rbind(bug.out.unm_conf_McCandless$sims.array[k,1,],bug.out.unm_conf_McCandless$sims.array[k,2,],bug.out.unm_conf_McCandless$sims.array[k,3,]))


# Summarizing results -----------------------------------------------------

# We can now summarized the MCMC results (median, percentiles 2.5th and 97.5th)
# medians and 95% credible intervals
res_unm_conf_McCandless <- round(
  t(
    apply(X=estimates_unm_conf_McCandless, MARGIN=2, FUN=quantile, probs=c(0.5, 0.025, 0.975))
  ), 
  digits=3)
res_unm_conf_McCandless

# Plot the posterior distributions for the adjusted OR using McCandless (black) vs Lash (blue) vs the biased OR (red)
plot(density(x=estimates_unm_conf_McCandless$ORquest), 
     main="Adjusted OR",
     xlab="Value", ylab="Density",
     ylim=c(0, 7),
     xlim=c(0, 0.7)
)

# Posterior distribution of OR_adj - Lash model
lines(density(estimates_unm_conf_lash$OR_adj), col="blue")

# Posterior distribution of OR_biased
lines(density(estimates_unm_conf_lash$OR_biased), col="red", lty="dashed")

