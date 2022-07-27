library(coda)
library(nimble)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

## Set seed
mySeed <- 0

#******************************************************************************#
#* DATA PREPARATION 
#******************************************************************************#

#-------------------------------------------------------------------------------
# !! IMPORTANT NOTE !! 
# The processed data required to run the model ('WildBoarIPM_Data.RData') is not
# openly available at the moment. 
# To request access to the data, please contact Dr. Marlène Gamelon: 
# marlene.gamelon@univ-lyon1.fr
#-------------------------------------------------------------------------------

## Load data
load('WildBoarIPM_Data.RData')

#  mydata_ter = Multi-state capture histories
#  Nb_Jeunes = Number of fetuses per pregnant female, by size class (harvested)
#  Nb_Gestantes = Number of pregnant females, by size class (harvested)
#  Nb_Part_Repro = Number of reproducing females (oestrus, ovulated, or pregnant) during hunting season, by size class (harvested)
#  Nb_Females = Number of harvested females for which reproductive status has been assessed


## Write size-at-harvest data
#  Numbers of females (small, medium, large) harvested between 1991 and 2016
yHS <- c(166, 94, 49, 85, 124, 198, 139, 116, 112, 177, 143, 170, 224, 173, 184, 186, 179, 265, 196, 164, 123, 177, 181, 119,  206, 148)
yHM <-c(56, 71, 31, 31, 75, 127, 53, 108, 100, 58, 68, 125, 99, 66, 78, 68, 187, 118, 64, 83, 108, 91, 31, 102, 69, 45)
yHL <- c(34, 34, 15, 20, 57, 90, 63, 51, 65, 39, 51, 46, 75, 31, 19, 36, 42, 38, 59, 42, 33, 32, 34, 31, 70, 33)


## Assemble Size-at-Harvest matrix
SaH <- rbind(yHS, yHM, yHL)
#dimnames(SaH) <- list(c("Small", "Medium", "Large"), c(1991:2016))


## Set number of years in size-at-harvest data
nyears <- length(yHS)


## Format reproduction data

# Size-class specific fetus numbers
nFetus = t(as.matrix(Nb_Jeunes)) # Number of foetus counted per year (rows) for small, medium and large females (columns)

# Size-class specific number pregnant
nPreg = t(as.matrix(Nb_Gestantes)) # Number of small, medium and large females (in columns) pregnant annually (in rows) 

# Size-class specific number reproducting
nRep = t(as.matrix(Nb_Part_Repro)) # Number of small, medium and large females (in columns) in oestrus, ovulated, or pregnant annually (in rows) 

# Size-class specific number sampled (for reproductive status)
nFemS = t(as.matrix(Nb_Females)) # Number of small, medium and large females (in columns) for which reproductive status has been assessed (in rows) 


## Format mark-recapture-recovery data
mydata <- unname(as.matrix(mydata_ter))
str(mydata)
table(mydata)

# Remove individual 39
mydata <- mydata[-39,]
# --> Gets dropped because it's harvested twice (error)

# Re-define states to reflect implemented CMRR model structure
# 1 = captured alive, small
# 2 = captured alive, medium
# 3 = captured alive, large
# 4 = reported dead (harvest), small
# 5 = reported dead (harvest), medium
# 6 = reported dead (harvest), large
# 7 = not observed

mydata[which(mydata==0)] <- 7

# Set number of individuals and capture occasions (years)
n <- dim(mydata)[[1]] 
K <- dim(mydata)[[2]] 

# Determine individual times of first and last capture
first <- rep(NA, n)
last <- rep(NA, n)

for (i in 1:n){
  first[i] <- min(which(mydata[i,] < 7))
  
  if(any(mydata[i,]%in%c(4:6))){
  	
  	last[i] <- min(which(mydata[i,]%in%c(4:6)))
  }else{
  	last[i] <- K
  }
}


## Make categorical acorn covariate
# 1 = N (no acorn production)
# 2 = A (average acorn production) 
# 3 = H (high acorn production)

AcornData <- data.frame(
  Year = 1983:2016, 
	AcornLevel = c(1,2,1,1,2,1,3,1,1,2,1,3,1,1,1,3,1,1,3,1,2,1,1,2,2,1,1,1,1,1,1,2,1,2)
  )

AcornData$BoarYear <- AcornData$Year + 1
# --> Adding +1 to get the "boar year", which starts with reproduction in March, which in turn is affected by Acorn availability in the preceding season (Sep[t-1] to Mar[t])

AcornData <- subset(AcornData, BoarYear %in% c(1991:2016))


## Arrange into data and constants for NIMBLE
WB.data <- list(SaH = SaH, 
				        y = as.matrix(mydata),
                nFetus = nFetus, nPreg = nPreg, 
				        nRep = nRep, nFemS = nFemS) 

WB.constants <- list(Tmax = nyears, Z = dim(SaH)[[1]],
					           n.ind = n, first = first, last = last,
                     AcornCat = AcornData$AcornLevel) 


#*******************************************************************************#
#* MODEL CODE
#*******************************************************************************#

# --> Model Timing: Growth -> Natural mortality -> Harvest mortality

WB.code <- nimbleCode({
  
  #######################
  # 1) POPULATION MODEL #
  #######################
  
  #-------------------#
  # 1.1) Reproduction #
  #-------------------#
  
  # Parameters:
  # Off[z,t]: Number of female offspring produced by mothers in size class z in year t
  # N[z,t]: Number of females in size class z in year t
  # pB[z,t]: Breeding proportion (probability) of females in size class z in year t
  # nF[z,t]: Fetus number of females in size class z in year t
  # YOY[t]: Number of female offspring recruiting into the population as young-og-the year (YOY) in year t
  # s0[t]: New offspring's probability to survive first non-harvest season (year t)

  for(t in 1:Tmax){
    for(z in 1:Z){
      
      # Number of offpring produced by mothers in size class z
      Off[z,t] ~ dpois(marN[z,t]*pB[z,t]*0.5*nF[z,t])    
    }
    
    # Number of produced offspring recruiting into the population    
    YOY[t] ~ dbin(s0[t], sum(Off[1:Z,t]))

  }
  
  Off[1:Z, Tmax+1] <- 0
  
  #-------------#
  # 1.2) Growth #
  #-------------#
  
  for(t in 1:Tmax){

  	# Size class 1 - grow vs. not grow
  	marN1_plus[t] ~ dbin(gP[1,t], marN[1,t] + YOY[t]) # Size 1 growing to any size > 1
  	marN_g[1,1,t] <- (marN[1,t] + YOY[t]) - marN1_plus[t] # Size 1 remaining size 1
  	
  	# Size class 1 - grow to size class 2 vs. 3 (given growth)
  	marN_g[1,3,t] ~ dbin(gSL[t], marN1_plus[t]) # Growing size 1 entering size 3
  	marN_g[1,2,t] <- marN1_plus[t] - marN_g[1,3,t] # Growing size 1 entering size 2
  	
  	# Size class 2 - grow to size class 3 vs. not grow
  	marN_g[2,3,t] ~ dbin(gP[2,t], marN[2,t]) # Size 2 growing to size 3
  	marN_g[2,2,t] <-  marN[2,t] - marN_g[2,3,t] # Size 2 remaining size 2
  	marN_g[2,1,t] <- 0
  	
  	# Size class 3 - no growth
  	marN_g[3,3,t] <- marN[3,t] # Size 3 remaining size 3
  	marN_g[3,1:2,t] <- 0	
  }
  
  ## NOTE:
  # Ideally, growth (state transition) would be modelled using a multinomial likelihood.
  # However, NIMBLE's standard samplers do not work for latent multinomial models (i.e. with an un-fixed "size") as per this date (Feb 2021).
  # We therefore use sequential binomial likelihoods as a work-around. 
  
  
  #-----------------------------------#
  # 1.3) Non-harvest season (Mar-Sep) #
  #-----------------------------------#
  
  for(t in 1:Tmax){
  	for(z in 1:Z){
  		octN_g[z,t] ~ dbin(sN[z,t], sum(marN_g[1:Z,z,t])) # Survivors
  	}
  }
  
  
  #-------------------------------#
  # 1.4) Harvest season (Oct-Mar) #
  #-------------------------------#
  
  for(t in 1:Tmax){
  	for(z in 1:Z){
  		marN[z,t+1] ~ dbin(sH[z,t], octN_g[z,t]) # Survivors
  		H[z,t+1] <- octN_g[z,t] - marN[z,t+1] # Harvests
  	}
  }
  
  H[1:Z,1] <- 0
  
  
  
  #############################
  # 2) SIZE-AT-HARVEST MODULE #
  #############################
  
  # Data:
  # SaH[z,t] = Size-at-Harvest matrix (number of size class z individuals reported as harvested in the Oct[t-1] to Feb[t] harvest season)
  
  # Parameters:
  # H[z,t] = number of size class z individuals harvested in the Oct[t-1] to Feb[t] harvest season
  # ll[z,t] = recovery rate of size class z individuals harvested in the Oct[t-1] to Feb[t] harvest season
  
  for(t in 2:Tmax){
    for(z in 1:Z){
      SaH[z,t] ~ dbin(ll[z,t], H[z,t])
    }
  }
  
  
  
  #####################################
  # 3) MARK-RECAPTURE-RECOVERY MODULE #
  #####################################
  
  # Data:
  # y[ind,t] = capture history of indivdual 'ind' (in each year t)
  
  # Parameters:
  # ps[i,j,t] = transition probability from state i to state j in year t
  # po[i,obs,t] = observation probability of state i as obs in year t
  
  #------------------------------#
  # 3.1) State transition matrix #
  #------------------------------#
  
  # States:
  # 1 = alive, small
  # 2 = alive, medium
  # 3 = alive, large
  # 4 = newly dead (harvest), small
  # 5 = newly dead (harvest), medium
  # 6 = newly dead (harvest), large
  # 7 = dead (includes newly dead from non-harvest causes)
  
  for(t in 1:Tmax){
    
    # Alive transitions (stochastic)
    
    ps[1,1,t] <- (1-gP[1,t])*sN[1,t]*sH[1,t] 
    ps[1,2,t] <- gP[1,t]*(1-gSL[t])*sN[2,t]*sH[2,t]
    ps[1,3,t] <- gP[1,t]*gSL[t]*sN[3,t]*sH[3,t]
    ps[1,4,t] <- (1-gP[1,t])*sN[1,t]*(1-sH[1,t])
    ps[1,5,t] <- gP[1,t]*(1-gSL[t])*sN[2,t]*(1-sH[2,t]) 
    ps[1,6,t] <- gP[1,t]*gSL[t]*sN[3,t]*(1-sH[3,t]) 
    ps[1,7,t] <- 1-sum(ps[1,1:6,t])
    
    ps[2,1,t] <- 0 
    ps[2,2,t] <- (1-gP[2,t])*sN[2,t]*sH[2,t]
    ps[2,3,t] <- gP[2,t]*sN[3,t]*sH[3,t] 
    ps[2,4,t] <- 0
    ps[2,5,t] <- (1-gP[2,t])*sN[2,t]*(1-sH[2,t]) 
    ps[2,6,t] <- gP[2,t]*sN[3,t]*(1-sH[3,t]) 
    ps[2,7,t] <- 1-sum(ps[2,1:6,t])
    
    ps[3,1,t] <- 0 
    ps[3,2,t] <- 0
    ps[3,3,t] <- sN[3,t]*sH[3,t]
    ps[3,4,t] <- 0 
    ps[3,5,t] <- 0
    ps[3,6,t] <- sN[3,t]*(1-sH[3,t])
    ps[3,7,t] <- 1-sum(ps[3,1:6,t])
    
    
    # Dead transitions (deterministic)
    
    ps[4,1:6,t] <- 0
    ps[4,7,t] <- 1
    
    ps[5,1:6,t] <- 0
    ps[5,7,t] <- 1
    
    ps[6,1:6,t] <- 0
    ps[6,7,t] <- 1
    
    ps[7,1:6,t] <- 0
    ps[7,7,t] <- 1
    
  }
  
  
  #------------------------#
  # 3.2) Obervation matrix #
  #------------------------#
  
  # Observations
  # 1 = captured alive, small
  # 2 = captured alive, medium
  # 3 = captured alive, large
  # 4 = reported dead (harvest), small
  # 5 = reported dead (harvest), medium
  # 6 = reported dead (harvest), large
  # 7 = not observed
  
  for(t in 2:Tmax){
    
    po[1,1,t] <- pp[1,t]
    po[1,2:6,t] <- 0
    po[1,7,t] <- 1-pp[1,t]
    
    po[2,1,t] <- 0
    po[2,2,t] <- pp[2,t]
    po[2,3:6,t] <- 0
    po[2,7,t] <- 1-pp[2,t]
    
    po[3,1:2,t] <- 0
    po[3,3,t] <- pp[3,t]
    po[3,4:6,t] <- 0
    po[3,7,t] <- 1-pp[3,t]
    
    po[4,1:3,t] <- 0
    po[4,4,t] <- ll[1,t]
    po[4,5:6,t] <- 0
    po[4,7,t] <- 1-ll[1,t]
    
    po[5,1:4,t] <- 0
    po[5,5,t] <- ll[2,t]
    po[5,6,t] <- 0
    po[5,7,t] <- 1-ll[2,t]
    
    po[6,1:5,t] <- 0
    po[6,6,t] <- ll[3,t]
    po[6,7,t] <- 1-ll[3,t]
    
    po[7,1:6,t] <- 0
    po[7,7,t] <- 1
    
  }
  
  
  #------------------------#
  # 3.3) Likelihood (dcat) #
  #------------------------#
  
    # Likelihood
  for(i in 1:n.ind){
  	
  	# Define latent state at first capture
  	x[i, first[i]] <- y[i, first[i]]
  	
  	for(t in (first[i]+1):last[i]){
  		
  		# State process: draw x(t) given x(t-1)
  		x[i, t] ~ dcat(ps[x[i, t-1], 1:7, t-1])
  		
  		# Observation process: draw y(t) given x(t)
  		y[i, t] ~ dcat(po[x[i, t], 1:7, t])
  	}
  	
  }
  
  
  
  #######################
  # 4) FECUNDITY MODULE #
  #######################
  
  # Data: 
  # nRep[z,t] = total number of size class z harvested females reproducing in year t
  # nFemS[z,t] = total number of size class z harvested females sampled for reproductive status in year t
  # nPreg[z,t] = total number of size class z harvested females pregnant in year t
  # nFetus[z,t] = total number of fetuses counted for pregnant, harvested size class z females in year t
  
  # Parameters:
  # pB[z,t] = breeding probability of size class z females in year t
  # nF[z,t] = fetus number of size class z females in year t
  
  for(t in 1:Tmax){
    for(z in 1:Z){
      
      #--------------------------------------#
      # 4.1) Breeding proportion/probability #
      #--------------------------------------#
      
      nRep[z,t] ~ dbin(pB[z,t], nFemS[z,t])	
      
           
      #-------------------#
      # 4.2) Fetus number #
      #-------------------#
      
      nFetus[z,t] ~ dpois(nPreg[z,t]*nF[z,t])
      
    }
  }
  
  
  
  #############################
  # 5) PRIORS AND CONSTRAINTS #
  #############################
  
  #---------------------------#
  # 5.1) Mortality parameters #
  #---------------------------#
  
  for(z in 1:Z){
    
    sN[z,1:Tmax] <- exp(-mN[z,1:Tmax])
    sH[z,1:Tmax] <- exp(-mH[z,1:Tmax])
    
    log(mN[z,1:Tmax]) <- log(Mu.mN[z]) + epsilon.mN[1:Tmax]
    log(mH[z,1:Tmax]) <- log(Mu.mH[z]) + epsilon.mH[1:Tmax]
    
    Mu.mN[z] ~ dunif(0, 5)
    Mu.mH[z] ~ dunif(0, 5)
    
  }
  
  for(t in 1:Tmax){
    
    epsilon.mN[t] ~ dnorm(0, sd = sigma.mN)
    epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
  }
  
  sigma.mN ~ dunif(0, 5)
  sigma.mH ~ dunif(0, 5)
  
  
  #------------------------#
  # 5.2) Growth parameters #
  #------------------------#
  
  for(z in 1:(Z-1)){
  	
  	logit(gP[z,1:Tmax]) <- logit(Mu.gP[z]) + epsilon.gP[z,1:Tmax]

  	Mu.gP[z] ~ dunif(0, 1)
  	sigma.gP[z] ~ dunif(0, 5)
  }
  
  
  for(t in 1:Tmax){
    
  	for(z in 1:(Z-1)){
  		epsilon.gP[z,t] ~ dnorm(0, sd = sigma.gP[z])
  	}
    
    epsilon.gSL[t] ~ dnorm(0, sd = sigma.gSL)
  }
  
  logit(gSL[1:Tmax]) <- logit(Mu.gSL) + epsilon.gSL[1:Tmax]
  
  Mu.gSL ~ dunif(0, 1)
  sigma.gSL ~ dunif(0, 5)
  
  
  #------------------------------------#
  # 5.3) Recapture/Recovery parameters #
  #------------------------------------#
  
  for(z in 1:Z){
    
    logit(pp[z,1:Tmax]) <- logit(Mu.pp[z]) + epsilon.pp[1:Tmax]
    
    logit(ll[z,1:Tmax]) <- logit(Mu.ll) + epsilon.ll[1:Tmax]
    
    Mu.pp[z] ~ dunif(0, 1)
  }
  
  Mu.ll ~ dunif(0, 1)

  
  for(t in 1:Tmax){
    
    epsilon.pp[t] ~ dnorm(0, sd = sigma.pp)
    epsilon.ll[t] ~ dnorm(0, sd = sigma.ll)
  }
  
  sigma.pp ~ dunif(0, 5)
  sigma.ll ~ dunif(0, 5)
  
  
  #------------------------------#
  # 5.3) Reproduction parameters #
  #------------------------------#
  
  for(z in 1:Z){
    
    for(t in 1:Tmax){
    	logit(pB[z,t]) <- logit(Mu.pB[z,AcornCat[t]]) + epsilon.pB[t]
    }
    
    for(a in 1:3){
   		Mu.pB[z,a] ~ dunif(0, 1)
    }

    log(nF[z,1:Tmax]) <- log(Mu.nF[z]) + epsilon.nF[1:Tmax]

    Mu.nF[z] ~ dunif(0, 10)
    
  }
  
  
  for(t in 1:Tmax){
    
    epsilon.pB[t] ~ dnorm(0, sd = sigma.pB)
    epsilon.nF[t] ~ dnorm(0, sd = sigma.nF)
  }
  
  sigma.pB ~ dunif(0, 5)
  sigma.nF ~ dunif(0, 5)
  
  
  #----------------------------#
  # 5.4) Population parameters #
  #----------------------------#
  
  ## Initial population sizes

  for(z in 1:Z){
  	initN[z] ~ T(dnorm(50, 50), 0, Inf)
  	marN[z,1] <- round(initN[z])
  }
  
  
  #----------------------------------------#
  # 5.5) Early survival ('Free' Parameter) #
  #----------------------------------------#
  
  for(t in 1:Tmax){
  	
  	s0[t] <- exp(-m0[t])
  	log(m0[t]) <- log(Mu.m0) + epsilon.m0[t]
   	epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
  }
  
  Mu.m0 ~ dunif(0, 5)
  sigma.m0 ~ dunif(0, 5)
  
})


#*******************************************************************************#
#* INITIAL VALUES
#*******************************************************************************#

## Load functions for simulating initial values
source('WildBoarIPM_InitialValuesSim.R')


#*******************************************************************************#
#* MCMC SETUP
#*******************************************************************************#

## Setting parameters to monitor
parameters <- c('Mu.mN', 'Mu.mH', 'sigma.mN', 'sigma.mH', 'mN', 'mH',
				'Mu.gP', 'Mu.gSL', 'gP', 'gSL',
				'Mu.pp', 'Mu.ll', 'sigma.pp', 'sigma.ll', 'pp', 'll',
				'Mu.nF', 'sigma.nF', 'Mu.pB', 'sigma.pB', 'pB', 'nF',
				'Mu.m0', 'sigma.m0', 'm0',
				'marN', 'marN_g', 'octN_g', 'H', 
				'Off', 'YOY'
				)

## MCMC settings
ni <- 150000
nb <- 70000
nt <- 10
nc <- 4

## Sample initial values (4 chains)
Inits <- list(WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5)), 
              WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5)), 
              WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5)),
              WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5)))


#*******************************************************************************#
#* RUN
#*******************************************************************************#

WB.IPM <- nimbleMCMC(code = WB.code, constants = WB.constants, data = WB.data, inits = Inits, monitors = parameters, niter = ni, nburnin = nb, nchains = nc, thin = nt, setSeed = mySeed, samplesAsCodaMCMC = TRUE)

#save(WB.IPM, file = 'WildBoarIPM_MCMCsamples.RData')
