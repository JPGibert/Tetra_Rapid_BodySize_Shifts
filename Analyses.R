#########################################################
####  			Jean-Philippe Gibert, PhD. 			 ####	
####	Code for Gibert, Han, Wieczynski & Yammine   ####		
#########################################################

## INIT Packages
library("dplyr")
library("mgcv")
library("FME")
library("multispatialCCM")
library("ggplot2")
library("ggridges")
library("RColorBrewer")
library("viridis")
library("lattice")
setwd("~/Desktop/JP/Papers_in_progress/JP_Tetra_Eco_Pheno_no_T")

## LOAD DATA
data <- read.csv("~/Desktop/JP/Papers_in_progress/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data.csv")

data_1 <- read.csv("~/Desktop/JP/Papers_in_progress/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_1.csv")
data_2 <- read.csv("~/Desktop/JP/Papers_in_progress/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_2.csv")
data_3 <- read.csv("~/Desktop/JP/Papers_in_progress/JP_Tetra_Eco_Pheno_no_T/Manuscript/Git_Hub/data_3.csv")

data <- rbind(data_1,data_2,data_3)

# Check data
head(data)

##-------------------------------------------------------------------------------
## 1) DENSITIES OVER TIME

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
## 2) TRAITS OVER TIME

# Summarize trait data by treatment, day and jar
trait_dat <- data %>%
  group_by(Treatment,Day, Jar) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE),SD = sd(Volume, na.rm=TRUE), MIN=min(Volume, na.rm=TRUE)) 

# Summarize trait data by treatment, day
trait_mean_dat <- data %>%
  group_by(Treatment, Day) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE))
 
# Plot!
plot(log(trait_dat$Vol) ~ trait_dat$Day,, pch=16)


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

# III) Data prepping for fitting	
	## This eliminates *one* outlier at Day 9 (well outside two SD from mean)
	## And also makes it so that t starts at t=0, as opposed to 1, which is needed for fitting.

# For Densities
dat_22S <- new_dat 	
dat_22S <- dat_22S[-which(dat_22S$Day==10 & dat_22S$Density<3000),] # Eliminates outlier at Day 9
dat_22S <- cbind(time = dat_22S$Day-1, R = dat_22S$Density, sd = rep(0.45,length(dat_22S$Density)))	

# For Traits 
dat_22Strait <- trait_dat
dat_22Strait <- cbind(time = dat_22Strait$Day-1, M = dat_22Strait$Vol, sd = rep(0.45,length(dat_22Strait$Vol)))	

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

## Run models with fitted parameters for plotting
finalnoK <- LV_model_noK(pars = c(coef(Fit_22SnoK), a=2.2*10^-4,b=8.22e+06)) 
finalSD2 <- LV_model_SD2(pars = c(coef(Fit_22SSD2), a=2.2*10^-4,b=8.22e+06))


##-------------------------------------------------------------------------------
## 3) PARAMETER UNCERTAINTY THROUGH MCMC 

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


##--------------------------------------------------------------------------------------------
## 4) MODEL SELECTION AND AIC

# NOTE: Because Model selection is done using the parameter values taken from the MCMC procedure, they can vary a little bit from those reported in th emain text. However, the main result (SD model better than EE model) has never been different in all my runs (hard to make up for a +40 AIC difference).

### I) SUPPLY AND DEMAND MODEL
## 1) Calculate the Sum or Squared Errors contribution
to_AIC <- LV_model_SD2(pars = c(gamma_mean, 
		eta_mean,m_mean,D_mean,i_mean, a=2.2*10^-4,b=8.22e+06))
to_AIC <- to_AIC %>%
	filter(time %in% seq(0,14,1)) %>%
	slice(rep(1:n(), each = 6))  # This lines repeats each row/time step 6 times

to_AIC_den <- to_AIC[-57,]# Gets rid of one line at day 9 to account for fact that one time point is missing from density time series
dens_AIC <- 89*log(sum((dat_22S[,2]-to_AIC_den[,2])^2))  #89 is the number of data points
trait_AIC <- 90*log(sum((dat_22Strait[,2]-to_AIC[,3])^2)) 

## 2) Calculate the AIC
AIC_SD <- 2*7 + dens_AIC + trait_AIC 

pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06)

### II) ECO-EVO MODEL
## 1) Calculate the Sum or Squared Errors contribution
to_AIC_EE <- LV_model_noK(pars = c(gamma=gammaEC_mean,eta=etaEC_mean,m=mEC_mean,sigma=sigmaEC_mean,a=2.2*10^-4,b=8.22e+06))
to_AIC_EE <- to_AIC_EE %>%
	filter(time %in% seq(0,14,1)) %>%
	slice(rep(1:n(), each = 6))  # This lines repeats each row/time step 6 times

to_AIC_den_EE <- to_AIC_EE[-57,]# Gets rid of one line at day 9 to account for fact that one time point is missing from density time series
dens_AIC_EE <- 89*log(sum((dat_22S[,2]-to_AIC_den_EE[,2])^2))  #89 is the number of data points
trait_AIC_EE <- 90*log(sum((dat_22Strait[,2]-to_AIC_EE[,3])^2)) 

## 2) Calculate the AIC
AIC_EE <- 2*6 + dens_AIC_EE + trait_AIC_EE 

## DeltaAIC
AIC_EE-AIC_SD


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

	## The right Embedding dimension is the smallest one that maximizes rho (from tutorials above)
## Based on rho vs E plot, we choose:
E_22Stable <- 5
E_TR22Stable <- 5

### FROM CCM MANUAL
#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
par(mfrow=c(2,2))
## 22STABLE
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_22Stable, Tr_TS_22Stable, E_22Stable, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_22Stable, Ab_TS_22Stable, E_TR22Stable, tau=1, iterations=800)
#Test for significant causal signal

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="22 Stable")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
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
	# Eliminate outlier at Day 9
dat_22S <- dat_22S[-which(dat_22S$Day==10 & dat_22S$Density<3000),] 
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
par(mfrow=c(3,2), mar=c(2,3,1,1))
#SD Model
for(i in 1:700){
	if(i==1){
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22S[,2]~dat_22S[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(2,10000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
		lines(finalSD2$time, finalSD2$R, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}else{ 
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		lines(to_plot_finalSD2$time, to_plot_finalSD2$R, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}		
}
lines(finalSD2$time, finalSD2$R, col=viridis(20)[8], lwd=3)
points(dat_22S[,2]~dat_22S[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))
			
for(i in 1:700){
	if(i==1){
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		lines(finalSD2$time, finalSD2$M, col=viridis(20)[8], lwd=2)
		lines(to_plot_finalSD2$time, to_plot_finalSD2$M, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}else{
		to_plot_finalSD2 <- LV_model_SD2(pars = c(gamma=rnorm(1,gamma_mean,gamma_sd), 
		eta=rnorm(1,eta_mean,eta_sd),m=rnorm(1,m_mean,0),D = rnorm(1,D_mean,D_sd),i=rnorm(1,i_mean,i_sd), a=2.2*10^-4,b=8.22e+06))	
		lines(to_plot_finalSD2$time, to_plot_finalSD2$M, col=viridis(20, alpha=0.07)[8], lwd=0.3)
	}		
}
lines(finalSD2$time, finalSD2$M, col=viridis(20)[8], lwd=3)
points(dat_22Strait[,2]~dat_22Strait[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))

# Eco-evo
for(i in 1:700){
	if(i==1){
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))	
		plot(dat_22S[,2]~dat_22S[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(2,10000),axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
		lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=2)
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}else{
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))		
		lines(to_plot_finalnoK$time, to_plot_finalnoK$R, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}		
}
lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=3)
points(dat_22S[,2]~dat_22S[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))
	
for(i in 1:700){
	if(i==1){
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))
		plot(dat_22Strait[,2]~dat_22Strait[,1], xlab = "", ylab = "", pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE)
		box(bty="l", lwd=2)
		axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
		axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
		lines(finalnoK$time, finalnoK$M, col=viridis(20)[15], lwd=2)
		lines(to_plot_finalnoK$time, to_plot_finalnoK$M, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}else{
		to_plot_finalnoK <- LV_model_noK(pars = c(gamma=rnorm(1,gammaEC_mean,gammaEC_sd), eta=rnorm(1,etaEC_mean,etaEC_sd),m=rnorm(1,mEC_mean,mEC_sd),sigma=rnorm(1,sigmaEC_mean,sigmaEC_sd), a=2.2*10^-4,b=8.22e+06))
		lines(to_plot_finalnoK$time, to_plot_finalnoK$M, col=viridis(20, alpha=0.07)[15], lwd=0.3)
	}		
}
lines(finalnoK$time, finalnoK$M, col=viridis(20)[15], lwd=3)
points(dat_22Strait[,2]~dat_22Strait[,1], xlab = "time", ylab = "N", pch=1, col="darkgrey", ylim=c(2,10000))


## Panels e & f

#par(mfrow=c(1,2), mar=c(2,3,1,1))
plot(dat_22S[,2]~dat_22S[,1], pch=1, col="darkgrey", ylim=c(2,10000), axes=FALSE, ylab="",xlab="")
box(bty="l", lwd=2)
axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
axis(2,c(0,2,4,6,8,10),at=seq(0,10000,2000), tck=0.010, lwd=2, hadj=0.2,  las=TRUE)
lines(finalSD2$time, finalSD2$R, col=viridis(20)[8], lwd=2)
lines(finalnoK$time, finalnoK$R, col=viridis(20)[15], lwd=2)
 
plot(dat_22Strait[,2]~dat_22Strait[,1], pch=1, col="darkgrey", ylim=c(4000,20000), axes=FALSE, ylab="",xlab="")
box(bty="l", lwd=2)
axis(1,seq(0,15,2), tck=0.010, lwd=2, padj=-1)
axis(2,c(5,10,15,20),at=seq(5000,20000,5000), tck=0.010, lwd=2, hadj=0.2, las=TRUE)
lines(finalSD2$time, finalSD2$M, col=viridis(20)[8], lwd=2)
lines(finalnoK$time, finalnoK$M, col=viridis(20)[15], lwd=2)


### THE END
