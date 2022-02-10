

## INIT Packages
library("dplyr")
library("mgcv")
library("FME")
library("multispatialCCM")
library("ggridges")
library("RColorBrewer")
library("viridis")
library("lattice")
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub")

## LOAD DATA

#1) Time Series Data
data_1 <- read.csv("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_1.csv")
data_2 <- read.csv("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_2.csv")
data_3 <- read.csv("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_3.csv")

data <- rbind(data_1,data_2,data_3)

#2) Experimental Manipulation Data
# Manipulation of Density Experimental Data
data_1 <- read.csv("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/x1_Functional_Ecology/2_Resubmission/New_analyses/Density_on_Size.csv")
# Manipulation of Body Size Experimental Data
data_2 <- read.csv("~/Desktop/JP/Papers_in_review_submitted/JP_Tetra_Eco_Pheno_no_T/Manuscript/x1_Functional_Ecology/2_Resubmission/New_analyses/Size_on_Density.csv")


#############################################################################################
## Time Series-related analyses (Figs 2, and 4)

##-------------------------------------------------------------------------------
## 1) CHANGE IN DENSITIES OVER TIME

# Summarize density data by treatment, day and jar
new_dat <- data %>%
  group_by(Treatment,Day, Jar) %>%
  summarize(Density = mean(Density, na.rm = TRUE))
print(new_dat, n=Inf)

# As before but only for means across Jars and coefficient of variation
means_new_dat <- data %>%
	group_by(Treatment, Day) %>%
	summarize(Abundance=mean(Density, na.rm=TRUE), CV=sd(Density, na.rm=TRUE)/mean(Density, na.rm=TRUE))

# Check data
plot(Density~jitter(Day,2), data=new_dat, pch=16)
	

##-------------------------------------------------------------------------------
## 2) CHANGE IN BODY SIZE OVER TIME

# Summarize trait data by treatment, day and jar
trait_dat <- data %>%
  group_by(Treatment,Day, Jar) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE),SD = sd(Volume, na.rm=TRUE), MIN=min(Volume, na.rm=TRUE)) 

# Summarize trait data by treatment, day
trait_mean_dat <- data %>%
  group_by(Treatment, Day) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE))
 
# Plot!
plot(trait_dat$Vol ~ trait_dat$Day,pch=16)


##-------------------------------------------------------------------------------
## 3) ODE FITTING

# I) ECO-EVO Model
# Define model to be fitted (Model with shifting body size)
LV_model_noK <- function (pars, R_0 = 10, M_0 = 14231) {
	derivs <- function(time, y, pars) {
		with (as.list(c(pars, y)), {
			dR <- a*(M^gamma)*R*(1-R/(b*(M^(eta))))-m*R
			dM <- sigma*a*(M^(gamma-eta-1))*(b*(M^eta)*gamma + (eta-gamma)*R)
			return(list(c(dR,dM),logR=log(R), logM=log(M)))
			}
		)
	}
# Initial conditions
y <- c(R = R_0, M = M_0)
times <- c(seq(0, 14, 0.1)) 	
out <- ode(y = y, parms = pars, times = times, func = derivs)
as.data.frame(out)
}

# II) Supp/Dem Model 
# Define model to be fitted (Model with shifting body size)
LV_model_SD2 <- function (pars, R_0 = 10, M_0 = 14231) {
	derivs <- function(time, y, pars) {
		with (as.list(c(pars, y)), {
			dR <- a*(M^gamma)*R*(1-R/(b*(M^(eta))))-m*R
			dM <- (i*(b*(M^(eta)))/R - D*M)
			return(list(c(dR,dM),logR=log(R), logM=log(M)))
			}
		)
	}
# Initial conditions
y <- c(R = R_0, M = M_0)
times <- c(seq(0, 14, 0.1)) 	
out <- ode(y = y, parms = pars, times = times, func = derivs)
as.data.frame(out)
}

# III) Eco-Evo + Supp/Dem Model 
# Define model to be fitted (Model with shifting body size)
LV_model_SDEE <- function (pars, R_0 = 10, M_0 = 14231) {
	derivs <- function(time, y, pars) {
		with (as.list(c(pars, y)), {
			dR <- a*(M^gamma)*R*(1-R/(b*(M^(eta))))-m*R
			dM <- (i*(b*(M^(eta)))/R - D*M) + (sigma*a*(M^(gamma-eta-1))*(b*(M^eta)*gamma + (eta-gamma)*R))
			return(list(c(dR,dM),logR=log(R), logM=log(M)))
			}
		)
	}
# Initial conditions
y <- c(R = R_0, M = M_0)
times <- c(seq(0, 14, 0.1)) 	
out <- ode(y = y, parms = pars, times = times, func = derivs)
as.data.frame(out)
}


# IV) Data prepping for fitting	
	## This makes it so that t starts at t=0, as opposed to 1, which is needed for fitting.

# For Densities
dat_22S <- new_dat 	
dat_22S <- cbind(time = dat_22S$Day-1, R = dat_22S$Density, sd = rep(0.45,length(dat_22S$Density)))	

# For Traits 
dat_22Strait <- trait_dat
dat_22Strait <- cbind(time = dat_22Strait$Day-1, M = dat_22Strait$Vol, sd = rep(0.45,length(dat_22Strait$Vol)))	

# For Biomass (Only important for figure S10 and associated discussion. Note: Biomass isn't fitted, it's just overlayed, fits are done with N and M). 
for_bio <- merge(trait_dat, new_dat, intersect=c("Treatment", "Day", "Jar"))
#Convert volume in microm^3 to cm^3 because density of water is 1g/cm^3, so as to estimate Biomass in g/mL
for_bio$Mass <- for_bio$Vol*10^-12 # Assuming density=1g/cm^3, then total volume equals total mass.
for_bio$Biomass <- for_bio$Density*for_bio$Mass
dat_22Sbiomass <- for_bio
dat_22Sbiomass <- cbind(time = dat_22Sbiomass$Day-1, M = dat_22Sbiomass$Biomass, sd = rep(0.45,length(dat_22Sbiomass$Biomass)))	


# IV) Model fitting using FME

	# 1) EcoEvo
# Define function wrapper with ODE model, datasets to be used, and cost function to be optimized
LVcost_22S_noK <- function (pars) {
	out <- LV_model_noK(pars)
	cost <- modCost(model = out, obs = dat_22S, err = "sd")
	return(modCost(model = out, obs = dat_22Strait, err = "sd", cost = cost))
}
#NOTE: Model parameters are already initilized close to their final values here so fitting converges
LVcostnoK <- function(Npars){LVcost_22S_noK(c(Npars, a=2.2*10^-4,b=8.22e+06))}
parsnoK <- c(gamma=0.907, eta=-0.26, sigma=0.809,m=1.08) 

## Actual model fitting using function modFit()
Fit_22SnoK <- modFit(f = LVcostnoK, p = parsnoK)

## Assess model fit and parameter p-values using FME's summary, then check values
summary(Fit_22SnoK)
coef(Fit_22SnoK)

	# 2) Supply-Demand
LVcost_22S_SD2 <- function (pars) {
	out <- LV_model_SD2(pars)
	cost <- modCost(model = out, obs = dat_22S, err = "sd")
	return(modCost(model = out, obs = dat_22Strait, err = "sd", cost = cost))
}
#NOTE: Model parameters are already initilized close to their final values here so fitting converges
LVcostSD2 <- function(Npars){LVcost_22S_SD2(c(Npars, a=2.2*10^-4,b=8.22e+06))}
parsSD2 <- c(gamma=0.9, eta=-0.23, m=0.2, D=0.08, i=24)

## Actual model fitting using function modFit()
Fit_22SSD2 <- modFit(f = LVcostSD2, p = parsSD2)

## Assess model fit and parameter p-values using FME's summary, then check values
summary(Fit_22SSD2)
coef(Fit_22SSD2)

	# 3) Mixed plasticity and eco-evo model
LVcost_22S_SDEE <- function (pars) {
	out <- LV_model_SDEE(pars)
	cost <- modCost(model = out, obs = dat_22S, err = "sd")
	return(modCost(model = out, obs = dat_22Strait, err = "sd", cost = cost))
}
#NOTE: Model parameters are already initilized close to their final values here so fitting converges
LVcostSDEE <- function(Npars){LVcost_22S_SDEE(c(Npars, a=2.2*10^-4,b=8.22e+06))}
parsSDEE <- c(gamma=0.9, eta=-0.23, m=0.87, D=0.08, i=1, sigma=1)

## Actual model fitting using function modFit()
Fit_22SSDEE <- modFit(f = LVcostSDEE, p = parsSDEE)

## Assess model fit and parameter p-values using FME's summary, then check values
summary(Fit_22SSDEE)
coef(Fit_22SSDEE) ## Sigma not significant!

## Run models with fitted parameters for plotting
finalnoK <- LV_model_noK(pars = c(coef(Fit_22SnoK), a=2.2*10^-4,b=8.22e+06)) 
finalSD2 <- LV_model_SD2(pars = c(coef(Fit_22SSD2), a=2.2*10^-4,b=8.22e+06))
finalSSDEE <- LV_model_SDEE(pars = c(coef(Fit_22SSDEE), a=2.2*10^-4,b=8.22e+06))

##-------------------------------------------------------------------------------
## 3) PARAMETER UNCERTAINTY VIA MARKOV-CHAIN MONTE CARLO 

## 1) Parameter estimtion via MCMC for ECO-EVO model
var0 <- Fit_22SnoK$var_ms_unweighted # Initial model variance set as the residual mean square returned by the summary function
cov0 <- summary(Fit_22SnoK)$cov.scaled  #The -scaled- parameter covariances returned from the summary function are used as estimate of the proposal covariances between parameters. The scaling itself is a bit obscure, and is typically referred to as being scaled "appropriately". But what is meant by that seems to be reliant on how fast we want the MCMC algorithm to converge. More can be found here;Gelman A, Varlin JB, Stern HS, Rubin DB (2004). Bayesian Data Analysis. 

MCMC_ecoevo <- modMCMC(f = LVcostnoK , p = Fit_22SnoK$par, niter = 5000, jump = cov0,var0 = var0, wvar0 = 0.1, updatecov = 50)

	# To see that MCMC algorithm converged
plot(MCMC_ecoevo, Full = TRUE)
	# To see correlation between parameters
pairs(MCMC_ecoevo, nsample = 1000) 

# gamma
gammaEC_mean <- mean(MCMC_ecoevo$pars[,1])
gammaEC_sd <- sd(MCMC_ecoevo$pars[,1])
#eta
etaEC_mean <- mean(MCMC_ecoevo$pars[,2])
etaEC_sd <- sd(MCMC_ecoevo$pars[,2])
#sigma
sigmaEC_mean <- mean(MCMC_ecoevo$pars[,3])
sigmaEC_sd <- sd(MCMC_ecoevo$pars[,3])
#m
mEC_mean <- mean(MCMC_ecoevo$pars[,4])
mEC_sd <- sd(MCMC_ecoevo$pars[,4])

finalnoK <- LV_model_noK(pars = c(MCMC_ecoevo$bestpar, a=2.2*10^-4,b=8.22e+06))

# Credible intervals
quantile(MCMC_ecoevo$pars[,1], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_ecoevo$pars[,2], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_ecoevo$pars[,3], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_ecoevo$pars[,4], probs = c(0.025,0.975), na.rm = FALSE)

## 2) Parameter uncertainty via MCMC for SD model
var0 <- Fit_22SSD2$var_ms_unweighted # Initial model variance set as the residual mean square returned by the summary function
cov0 <- summary(Fit_22SSD2)$cov.scaled  # Done as above 

MCMC <- modMCMC(f = LVcostSD2 , p = Fit_22SSD2$par, niter = 5000, jump = cov0,var0 = var0, wvar0 = 0.1, updatecov = 50)

	# To see that MCMC algorithm converged
plot(MCMC, Full = TRUE)
	# To see correlation between parameters
pairs(MCMC, nsample = 1000) # Spurious correlations are pretty common, so user beware!

# gamma
gamma_mean <- MCMC$bestpar[1]
gamma_sd <- sd(MCMC$pars[,1])
#eta
eta_mean <- MCMC$bestpar[2]
eta_sd <- sd(MCMC$pars[,2])
#m
m_mean <- MCMC$bestpar[3]
m_sd <- sd(MCMC$pars[,3])
#D
D_mean <- MCMC$bestpar[4]
D_sd <- sd(MCMC$pars[,4])
#i
i_mean <- MCMC$bestpar[5]
i_sd <- sd(MCMC$pars[,5])

finalSD2 <- LV_model_SD2(pars = c(MCMC$bestpar, a=2.2*10^-4,b=8.22e+06))

# Credible intervals
quantile(MCMC$pars[,1], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC$pars[,2], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC$pars[,3], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC$pars[,4], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC$pars[,5], probs = c(0.025,0.975), na.rm = FALSE)


## 3) Parameter uncertainty via MCMC for SDEE model
var0 <- Fit_22SSDEE$var_ms_unweighted # Initial model variance set as the residual mean square returned by the summary function
cov0 <- summary(Fit_22SSDEE)$cov.scaled  # Done as above 

MCMC_SDEE <- modMCMC(f = LVcostSDEE , p = Fit_22SSDEE$par, niter = 5000, jump = cov0,var0 = var0, wvar0 = 0.1, updatecov = 50)

	# To see that MCMC algorithm converged
plot(MCMC_SDEE, Full = TRUE)
	# To see correlation between parameters
pairs(MCMC_SDEE, nsample = 1000) # Spurious correlations are pretty common, so user beware!

# gamma
gamma_SDEE_mean <- MCMC_SDEE$bestpar[1]
gamma_SDEE_sd <- sd(MCMC_SDEE$pars[,1])
#eta
eta_SDEE_mean <- MCMC_SDEE$bestpar[2]
eta_SDEE_sd <- sd(MCMC_SDEE$pars[,2])
#m
m_SDEE_mean <- MCMC_SDEE$bestpar[3]
m_SDEE_sd <- sd(MCMC_SDEE$pars[,3])
#D
D_SDEE_mean <- MCMC_SDEE$bestpar[4]
D_SDEE_sd <- sd(MCMC_SDEE$pars[,4])
#i
i_SDEE_mean <- MCMC_SDEE$bestpar[5]
i_SDEE_sd <- sd(MCMC_SDEE$pars[,5])
#sigma
sigma_SDEE_mean <- mean(MCMC_SDEE$pars[,3])
sigma_SDEE_sd <- sd(MCMC_SDEE$pars[,3])

finalSDEE <- LV_model_SDEE(pars = c(MCMC_SDEE$bestpar, a=2.2*10^-4,b=8.22e+06))

# Credible intervals
quantile(MCMC_SDEE$pars[,1], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_SDEE$pars[,2], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_SDEE$pars[,3], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_SDEE$pars[,4], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_SDEE$pars[,5], probs = c(0.025,0.975), na.rm = FALSE)
quantile(MCMC_SDEE$pars[,6], probs = c(0.025,0.975), na.rm = FALSE)

##--------------------------------------------------------------------------------------------
## 4) MODEL SELECTION AND AIC

# NOTE: Because Model selection is done using the parameter values taken from the MCMC procedure, they can vary a little bit from those reported in the main text. However, the main result (SD model better than EE model) has never been different in all my runs (hard to make up for a +40 AIC difference).

### I) SUPPLY AND DEMAND MODEL
## 1) Calculate the Sum or Squared Errors contribution
to_AIC <- LV_model_SD2(pars = c(gamma_mean, 
		eta_mean,m_mean,D_mean,i_mean, a=2.2*10^-4,b=8.22e+06))
to_AIC <- to_AIC %>%
	filter(time %in% seq(0,14,1)) %>%
	slice(rep(1:n(), each = 6))  # This lines repeats each row/time step 6 times

to_AIC_den <- to_AIC
dens_AIC <- 90*log(sum((dat_22S[,2]-to_AIC_den[,2])^2))  #89 is the number of data points
trait_AIC <- 90*log(sum((dat_22Strait[,2]-to_AIC[,3])^2)) 

## 2) Calculate the AIC
AIC_SD <- 2*7 + dens_AIC + trait_AIC ## As 2*k + AIC of density and trait model iwth k=number of parameters.


### II) ECO-EVO MODEL
## 1) Calculate the Sum or Squared Errors contribution
to_AIC_EE <- LV_model_noK(pars = c(gamma=gammaEC_mean,eta=etaEC_mean,m=mEC_mean,sigma=sigmaEC_mean,a=2.2*10^-4,b=8.22e+06))
to_AIC_EE <- to_AIC_EE %>%
	filter(time %in% seq(0,14,1)) %>%
	slice(rep(1:n(), each = 6))  # This lines repeats each row/time step 6 times

to_AIC_den_EE <- to_AIC_EE
dens_AIC_EE <- 90*log(sum((dat_22S[,2]-to_AIC_den_EE[,2])^2))  #89 is the number of data points
trait_AIC_EE <- 90*log(sum((dat_22Strait[,2]-to_AIC_EE[,3])^2)) 

## 2) Calculate the AIC
AIC_EE <- 2*6 + dens_AIC_EE + trait_AIC_EE 


### III) ECOEVO + SUPPLY AND DEMAND MODEL
## 1) Calculate the Sum or Squared Errors contribution
to_AIC_SDEE <- LV_model_SDEE(pars = c(gamma_SDEE_mean, 
		eta_SDEE_mean,m_SDEE_mean,D_SDEE_mean,i_SDEE_mean,sigma=sigma_SDEE_mean, a=2.2*10^-4,b=8.22e+06))
to_AIC_SDEE <- to_AIC_SDEE %>%
	filter(time %in% seq(0,14,1)) %>%
	slice(rep(1:n(), each = 6))  # This lines repeats each row/time step 6 times

to_AIC_den_SDEE <- to_AIC_SDEE
dens_AIC_SDEE <- 90*log(sum((dat_22S[,2]-to_AIC_den_SDEE[,2])^2))  #89 is the number of data points
trait_AIC_SDEE <- 90*log(sum((dat_22Strait[,2]-to_AIC_SDEE[,3])^2)) 

## 2) Calculate the AIC
AIC_SD_SDEE <- 2*8 + dens_AIC_SDEE + trait_AIC_SDEE ## As 2*k + AIC of density and trait model iwth k=number of parameters.

## DeltaAIC
AIC_EE-AIC_SD
AIC_SD_SDEE-AIC_SD

## AIC weights
exp(-0.5*(AIC_SD-AIC_SD))/(1+exp(-0.5*(AIC_EE-AIC_SD))+exp(-0.5*(AIC_SD_SDEE-AIC_SD)))
exp(-0.5*(AIC_EE-AIC_SD))/(1+exp(-0.5*(AIC_EE-AIC_SD))+exp(-0.5*(AIC_SD_SDEE-AIC_SD)))
exp(-0.5*(AIC_SD_SDEE-AIC_SD))/(1+exp(-0.5*(AIC_EE-AIC_SD))+exp(-0.5*(AIC_SD_SDEE-AIC_SD)))


##--------------------------------------------------------------------------------------------
## 5) CONVERGENT CROSS MAPPING

## I) PREP data for CCM

# Make mock time series for package
Ab_TS_22Stable <- rep(0,15*6+6)
Tr_TS_22Stable <- rep(0,15*6+6)

for(i in 1:6){
	# DENSITY DATA
	Dat_22Stable <- new_dat %>% 
	filter(Treatment=="22_Stable" & Jar==i)

	Ab_TS_22Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dat_22Stable$Density
	
	if(i == 6){
		Ab_TS_22Stable <- Ab_TS_22Stable[-1]
	}
	
	# TRAIT DATA
	Dattr_22Stable <- trait_dat %>% 
	filter(Treatment=="22_Stable" & Jar==i)
	Tr_TS_22Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dattr_22Stable$Vol
	if(i == 6){
		Tr_TS_22Stable <- Tr_TS_22Stable[-1]
	}
}

# Each replicated time series needs to be separated by an "NA" as per multispatialCCM package
	Ab_TS_22Stable[which(Ab_TS_22Stable==0)] <- NA
	Tr_TS_22Stable[which(Tr_TS_22Stable==0)] <- NA
	
## II) Find Embedding Dimension as per:
https://ha0ye.github.io/rEDM/articles/rEDM.html
https://mathbio.github.io/edmTutorials/ccm.html

#Calculate optimal E (Embedding Dimension)
maxE<-12 #Maximum E to test
#Matrix for storing output
Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("22Stable", "TR22Stable")
#Loop over potential E values and calculate predictive ability
#of each process for its own dynamics
for(E in 2:maxE) {
#Uses defaults of looking forward one prediction step (predstep)
#And using time lag intervals of one time step (tau)
	Emat[E-1,"22Stable"]<-SSR_pred_boot(A=Ab_TS_22Stable, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"TR22Stable"]<-SSR_pred_boot(A=Tr_TS_22Stable, E=E, predstep=1, tau=1)$rho
	}

#Look at plots to find E for each process at which
#predictive ability rho is maximized
matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("22Stable", "TR22Stable"), lty=1:2, col=1:2, lwd=2, bty="n") 
#Adding dashed line for clarity
abline(v=5,lty=2)

	## The right Embedding dimension is the smallest one that maximizes rho (from tutorials above)
## Based on rho vs E plot, we choose:
E_22Stable <- 5
E_TR22Stable <- 5

### FROM CCM MANUAL
#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)

## 22STABLE
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_22Stable, Tr_TS_22Stable, E_22Stable, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_22Stable, Ab_TS_22Stable, E_TR22Stable, tau=1, iterations=800)

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="22 Stable")

#Add +/- 1 standard deviation
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)

#Add +/- 1 standard deviation
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)
legend("topleft",
c("Density causes Traits", "Traits causes density"),
lty=c(1,2), col=c(1,2), lwd=2, bty="n")


##--------------------------------------------------------------------------------------------
## FIGURE 2

## 1) Data prepping for fitting
dat_22S <- new_dat %>%
	filter(Treatment == "22_Stable")

dat_22S <- cbind(time = dat_22S$Day-1, R = dat_22S$Density, sd = rep(0.45,length(dat_22S$Density)), Jar=as.numeric(dat_22S$Jar))	
# Using traits here: 
dat_22Strait <- trait_dat %>% 
	filter(Treatment == "22_Stable")
dat_22Strait <- cbind(time = dat_22Strait$Day-1, M = dat_22Strait$Vol, sd = rep(0.45,length(dat_22Strait$Vol)),Jar=as.numeric(dat_22Strait$Jar))	

dat_22S <- as.data.frame(dat_22S)
dat_22Strait <- as.data.frame(dat_22Strait)

dev.new()
## Actual Plotting
par(mfrow=c(1,3), mar=c(2,4,1,1))

# PANEL A
plot(dat_22S$R[which(dat_22S$Jar==1)]~dat_22S$time[which(dat_22S$Jar==1)], type="b", col="orange", pch=16, xlab = "", ylab = "", ylim=c(2,10000),axes=FALSE)
for(i in 2:6){points(dat_22S$R[which(dat_22S$Jar==i)]~dat_22S$time[which(dat_22S$Jar==i)], type="b", col="orange", pch=16)}
box(bty="l", lwd=2)
axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1, cex.axis=1.15)
axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE, cex.axis=1.15)

# PANEL B
plot(dat_22Strait$M[which(dat_22Strait$Jar==1)]~dat_22Strait$time[which(dat_22Strait$Jar==1)], type="b", col="purple", pch=16, xlab = "", ylab = "", ylim=c(4000,21000),axes=FALSE)
for(i in 2:6){points(dat_22Strait$M[which(dat_22Strait$Jar==i)]~dat_22Strait$time[which(dat_22Strait$Jar==i)], type="b", col="purple", pch=16)}
box(bty="l", lwd=2)
axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1, cex.axis=1.15)
axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE, cex.axis=1.15)

# PANEL C
#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col="orange", lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),xlab="", ylab="", axes=FALSE)
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col="orange")
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col="purple", lty=1, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col="purple")
box(lwd=2, bty='l')
axis(1,at=seq(0,80,20), tck=0.010, cex.axis=1.15, lwd.ticks=2,mgp=c(3, .5, 0))
axis(2,at=seq(0,1,0.2), tck=0.010, las=TRUE, cex.axis=1.15,lwd.ticks=2,mgp=c(3, .5, 0))


##--------------------------------------------------------------------------------------------
## FIGURE 3

dev.new()

par(mfrow=c(3,3), mar=c(2,3,1,1))
to_plot_finalSD2_mean <- LV_model_SD2(pars = c(gamma=gamma_mean[[1]],eta=eta_mean[[1]],m=m_mean[[1]],D = D_mean[[1]],i=i_mean[[1]], a=2.2*10^-4,b=8.22e+06))
## 1)
#SD MODEL
for(i in 1:700){
	if(i==1){
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22S[,2]~dat_22S[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(2,10000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
		#lines(finalSD2$time, finalSD2$R, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}else{ 
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}		
}
lines(to_plot_finalSD2_mean$time, to_plot_finalSD2_mean$R, col=viridis(20)[8], lwd=3)
points(dat_22S[,2]~dat_22S[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))
			
for(i in 1:700){
	if(i==1){
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalSD2$time, finalSD2$M, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalSD2$time, to_plot_finalSD2$M, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}else{
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		lines(to_plot_finalSD2$time, to_plot_finalSD2$M, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}		
}
lines(to_plot_finalSD2_mean$time, to_plot_finalSD2_mean$M, col=viridis(20)[8], lwd=3)
points(dat_22Strait[,2]~dat_22Strait[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))

for(i in 1:700){
	if(i==1){
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "", ylab = "", pch=1, col="darkgrey", axes=FALSE, ylim=c(10^-7,1.1*10^-4))
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(2,4,6,8,10),at=seq(2*10^-5,10^-4,2*10^-5), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalSD2$time, finalSD2$R*finalSD2$M*10^-12, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R*to_plot_finalSD2$M*10^-12, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}else{
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R*to_plot_finalSD2$M*10^-12, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}		
}
lines(to_plot_finalSD2_mean$time, to_plot_finalSD2_mean$R*to_plot_finalSD2_mean$M*10^-12, col=viridis(20)[8], lwd=3)
points(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))

## 2)
## ECO-EVO MODEL
to_plot_finalnoK_mean <- LV_model_noK(pars = c(gamma=gammaEC_mean[[1]], eta=etaEC_mean[[1]],m=mEC_mean[[1]],sigma=sigmaEC_mean[[1]], a=2.2*10^-4,b=8.22e+06))
for(i in 1:700){
	if(i==1){
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22S[,2]~dat_22S[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(2,10000),axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
		#lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=2)
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R, col=viridis(20, alpha=0.07)[12], lwd=0.3)
	}else{
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))		
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R, col=viridis(20, alpha=0.07)[12], lwd=0.3)
	}		
}
lines(to_plot_finalnoK_mean$time, to_plot_finalnoK_mean$R, col=viridis(20)[12], lwd=3)
points(dat_22S[,2]~dat_22S[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))
	
for(i in 1:700){
	if(i==1){
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))
		plot(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalnoK$time, finalnoK$M, col=viridis(20)[15], lwd=2)
		lines(to_plot_finalnoK$time, to_plot_finalnoK$M, col=viridis(20, alpha=0.07)[12], lwd=0.3)
	}else{
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))
		lines(to_plot_finalnoK$time, to_plot_finalnoK$M, col=viridis(20, alpha=0.07)[12], lwd=0.3)
	}		
}
lines(to_plot_finalnoK_mean$time, to_plot_finalnoK_mean$M, col=viridis(20)[12], lwd=3)
points(dat_22Strait[,2]~dat_22Strait[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))

for(i in 1:700){
	if(i==1){
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))		
		plot(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "", ylab = "", pch=1, col="darkgrey", axes=FALSE, ylim=c(10^-7,1.1*10^-4))
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(2,4,6,8,10),at=seq(2*10^-5,10^-4,2*10^-5), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalnoK$time, finalnoK$R*finalnoK$M*10^-12, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R*to_plot_finalnoK$M*10^-12, col=viridis(20, alpha=0.07)[12], lwd=0.3)	
	}else{
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R*to_plot_finalnoK$M*10^-12, col=viridis(20, alpha=0.07)[12], lwd=0.3)
	}		
}
lines(to_plot_finalnoK_mean$time, to_plot_finalnoK_mean$R*to_plot_finalnoK_mean$M*10^-12, col=viridis(20)[12], lwd=3)
points(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))


## 3)
## SD-ECO-EVO MODEL
to_plot_final_SDEE_mean <- LV_model_SDEE(pars = c(gamma=gamma_SDEE_mean[[1]], eta=eta_SDEE_mean[[1]],m=m_SDEE_mean[[1]],sigma=sigma_SDEE_mean[[1]],i=i_SDEE_mean[[1]],D=D_SDEE_mean[[1]], a=2.2*10^-4,b=8.22e+06))
for(i in 1:700){
	if(i==1){
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22S[,2]~dat_22S[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(2,10000),axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
		#lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=2)
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$R, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}else{
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))		
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$R, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}		
}
lines(to_plot_final_SDEE_mean$time, to_plot_final_SDEE_mean$R, col=viridis(20)[15], lwd=3)
points(dat_22S[,2]~dat_22S[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))

for(i in 1:700){
	if(i==1){
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=2)
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$M, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}else{
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))		
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$M, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}		
}
lines(to_plot_final_SDEE_mean$time, to_plot_final_SDEE_mean$M, col=viridis(20)[15], lwd=3)
points(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey")

for(i in 1:700){
	if(i==1){
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))		
		plot(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "", ylab = "", pch=1, col="darkgrey", axes=FALSE, ylim=c(10^-7,1.1*10^-4))
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(2,4,6,8,10),at=seq(2*10^-5,10^-4,2*10^-5), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		#lines(finalnoK$time, finalnoK$R*finalnoK$M*10^-12, col=viridis(20)[8], lwd=2)
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$R*to_plot_final_SDEE$M*10^-12, col=viridis(20, alpha=0.07)[8], lwd=0.3)	
	}else{
		to_plot_final_SDEE <- LV_model_SDEE(pars = c(gamma=rnorm(1,gamma_SDEE_mean,gamma_SDEE_sd), eta=rnorm(1,eta_SDEE_mean,eta_SDEE_sd),m=rnorm(1,m_SDEE_mean,m_SDEE_sd),sigma=rnorm(1,sigma_SDEE_mean,sigma_SDEE_sd),i=rnorm(1,i_SDEE_mean,i_SDEE_sd),D=rnorm(1,D_SDEE_mean,D_SDEE_sd), a=2.2*10^-4,b=8.22e+06))
		lines(to_plot_final_SDEE$time, to_plot_final_SDEE$R*to_plot_final_SDEE$M*10^-12, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}		
}
lines(to_plot_final_SDEE_mean$time, to_plot_final_SDEE_mean$R*to_plot_final_SDEE_mean$M*10^-12, col=viridis(20)[15], lwd=3)
points(dat_22Sbiomass[,2]~dat_22Sbiomass[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))



#############################################################################################
## Experimental Manipulation-related analyses (Figure 3)
library("effectsize")
## Density ==> Size
# Summarize trait data by treatment, day and jar
data_summ <- data_1 %>%
  group_by(Treatment, Day, Jar) %>%
  mutate(Day=factor(Day)) %>%
  summarize(Vol = mean(Biovolume..P..Spheroid., na.rm = TRUE), Density = mean(Ind.ml, na.rm = TRUE)) 

# Appendix Fig S2 
boxplot(Density~Treatment*Day, data=data_summ)
# Main text Fig 2a
boxplot(Vol~Treatment*Day, data=data_summ, range=3)

# Significant Interactive effect of Density and Day
summary(mod_1D <- lm(Vol~Treatment*Day, data=data_summ))

# Standardized interaction parameter
## method "refit" is the best for models with interactions (https://easystats.github.io/effectsize/articles/standardize_parameters.html) and two_sd=TRUE allow for direct comparison of continuous and discrete variable predictors (Gellman 2008).
standardize_parameters(mod_1D, method="refit", two_sd = TRUE)


## Body Size ==> Density
data_summ2 <- data_2 %>%
  group_by(Treatment, Day, Jar) %>%
   mutate(Day=factor(Day)) %>%
  summarize(Vol = mean(Biovolume..P..Spheroid., na.rm = TRUE), Density = mean(Ind.ml, na.rm = TRUE)) 

boxplot(Density~Treatment*Day, data=data_summ2, range=3)
boxplot(Vol~Treatment*Day, data=data_summ2)

summary(mod_1V <- lm(Vol~Treatment*Day, data=data_summ))
summary(mod_2V <- lm(Density~Treatment*Day, data=data_summ2))

standardize_parameters(mod_2V, method="refit", two_sd = TRUE)


## Calculate the magnitude of the imposed differences in SDs for Density and Volume
## For Density
data_init_dens <- data_summ %>%
	filter(Day==1)	
	
## For Volume
data_init_vol <- data_summ2 %>%
	filter(Day==1)	

data_init_vol$stdVol <- (data_init_vol$Vol-mean(data_init_vol$Vol))/sd(data_init_vol$Vol)
data_init_dens$stdDensity <- (data_init_dens$Density-mean(data_init_dens$Density))/sd(data_init_dens$Density)

# The magnitude of the difference will be given by
stdVol <-data_init_vol %>% 
	group_by(Treatment) %>%
	summarize(Vol=mean(stdVol))

stdDens <- data_init_dens %>% 
	group_by(Treatment) %>%
	summarize(Dens=mean(stdDensity))

## We then divide the interaction coefficient by that value to get an "Effect size divided by differnce in SD"
standardize_parameters(mod_1D, method="refit", two_sd = TRUE)[[2]][4]/(stdDens[1,2]-stdDens[2,2])
standardize_parameters(mod_2V, method="refit", two_sd = TRUE)[[2]][4]/(stdVol[1,2]-stdVol[2,2])
	# Pretty comparable suggesting strong effects of both.

## MAGNITUDE OF ALL DIFFERENCES USING TUKEY TEST
## 1) Main text results
## Body-size --> Density
## Test for significance of differences in Day 2 Using Tukey post-hoc
Dens.lm <- lm(Density~Treatment*Day, data=data_summ2)
Dens.av <- aov(Dens.lm)
summary(Dens.av)
tukey.test <- TukeyHSD(Dens.av, factor="Treatment") 
tukey.test # Check: Diff at Day 1 not sig

## Density --> Body-Size
## Test for significance of differences in Day 2 Using Tukey post-hoc
Vol.lm <- lm(Vol~Treatment*Day, data=data_summ)
Vol.av <- aov(Vol.lm)
summary(Vol.av)
tukey.test.Vol <- TukeyHSD(Vol.av, factor="Treatment")
tukey.test.Vol # Check: Diff at Day 1 not sig


## 2) Appendix Figs S2 and S3, Tukey test for differences at Day 0 (did treatments work?)
## Body-size --> Density
## Test for significance of differences in Day 2 Using Tukey post-hoc
Dens.lm_app <- lm(Vol~Treatment*Day, data=data_summ2)
Dens.av_app <- aov(Dens.lm_app)
summary(Dens.av_app)
tukey.test <- TukeyHSD(Dens.av_app, factor="Treatment") 
tukey.test # Check: Diff at Day 1 not sig

## Density --> Body-Size
## Test for significance of differences in Day 2 Using Tukey post-hoc
Vol.lm_app <- lm(Density~Treatment*Day, data=data_summ)
Vol.av_app <- aov(Vol.lm_app)
summary(Vol.av_app)
tukey.test.Vol <- TukeyHSD(Vol.av_app, factor="Treatment")
tukey.test.Vol # Check: Diff at Day 1 not sig





### THE END



