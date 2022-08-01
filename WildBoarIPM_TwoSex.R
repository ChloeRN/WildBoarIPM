library(coda)
library(nimble)
nimbleOptions(disallow_multivariate_argument_expressions = FALSE)

## Set seed
mySeed <- 20

## Set switch for model testing (TRUE = model run with only a subset CMRR data)
modelTest <- TRUE
#modelTest <- FALSE

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

# Numbers of males (small, medium, large) harvested between 1991 and 2016)
# NOTE: Data assumed identical to data for females for now
yHS_m <- yHS
yHM_m <- yHM
yHL_m <- yHL


## Assemble two-sex Size-at-Harvest matrix
SaH <- array(NA, dim = c(3, length(yHS), 2))
SaH[,,1] <- rbind(yHS, yHM, yHL)
SaH[,,2] <- rbind(yHS_m, yHM_m, yHL_m)
dimnames(SaH) <- list(c("Small", "Medium", "Large"), c(1991:2016), c('F', 'M'))


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

# Trim CMRR data (for model testing only)
if(modelTest){
  select.ind <- sample(1:dim(mydata)[1], 500, replace = F)
  mydata <- mydata[select.ind,]
}

# Duplicate CMRR data to simulate (identical) male capture histories
mydata <- rbind(mydata, mydata)

# Make vector for sex information for CMRR data
CH.sex <- rep(c(1, 2), each = dim(mydata)[1]/2)

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
					           n.ind = n, first = first, last = last, CH.sex = CH.sex,
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

  ## Linked model
  for(t in 1:Tmax){
    for(z in 1:Z){

      # Number of offpring (both sexes) produced by mothers in size class z
      Off_tot[z,t] ~ dpois(marN[z,t,1]*pB[z,t]*nF[z,t])
    }

    # Number of offspring of each sex
    Off[t,1] ~ dbin(0.5, sum(Off_tot[1:Z,t])) # Females
    Off[t,2] <- sum(Off_tot[1:Z,t]) - Off[t,1] # Males

    # Number of produced offspring recruiting into the population
    for(s in 1:2){
      YOY[t,s] ~ dbin(s0[t,s], Off[t,s])
    }
  }
  
  ## Unlinked model (treats m and f as two independent f populations)
  # for(s in 1:2){
  #   for(t in 1:Tmax){
  #     for(z in 1:Z){
  #       
  #       # Number of offpring produced by parents in size class z
  #       Off_tot[z,t,s] ~ dpois(marN[z,t,s]*pB[z,t]*nF[z,t]*0.5)    
  #     }
  #     
  #     # Number of produced offspring recruiting into the population
  #     YOY[t,s] ~ dbin(s0[t,s], sum(Off_tot[1:Z,t,s]))
  #   }
  # }

  
  #-------------#
  # 1.2) Growth #
  #-------------#
  
  for(s in 1:2){
    for(t in 1:Tmax){
      
      # Size class 1 - grow vs. not grow
      marN1_plus[t,s] ~ dbin(gP[1,t,s], marN[1,t,s] + YOY[t,s]) # Size 1 growing to any size > 1
      marN_g[1,1,t,s] <- (marN[1,t,s] + YOY[t,s]) - marN1_plus[t,s] # Size 1 remaining size 1
      
      # Size class 1 - grow to size class 2 vs. 3 (given growth)
      marN_g[1,3,t,s] ~ dbin(gSL[t,s], marN1_plus[t,s]) # Growing size 1 entering size 3
      marN_g[1,2,t,s] <- marN1_plus[t,s] - marN_g[1,3,t,s] # Growing size 1 entering size 2
      
      # Size class 2 - grow to size class 3 vs. not grow
      marN_g[2,3,t,s] ~ dbin(gP[2,t,s], marN[2,t,s]) # Size 2 growing to size 3
      marN_g[2,2,t,s] <-  marN[2,t,s] - marN_g[2,3,t,s] # Size 2 remaining size 2
      marN_g[2,1,t,s] <- 0
      
      # Size class 3 - no growth
      marN_g[3,3,t,s] <- marN[3,t,s] # Size 3 remaining size 3
      #marN_g[3,1:2,t,s] <- 0
      marN_g[3,2,t,s] <- 0
      marN_g[3,1,t,s] <- 0
    }
  }
 
  
  ## NOTE:
  # Ideally, growth (state transition) would be modelled using a multinomial likelihood.
  # However, NIMBLE's standard samplers do not work for latent multinomial models (i.e. with an un-fixed "size") as per this date (Feb 2021).
  # We therefore use sequential binomial likelihoods as a work-around. 
  
  
  #-----------------------------------#
  # 1.3) Non-harvest season (Mar-Sep) #
  #-----------------------------------#
  
  for(s in 1:2){
    for(t in 1:Tmax){
      for(z in 1:Z){
        octN_g[z,t,s] ~ dbin(sN[z,t,s], sum(marN_g[1:Z,z,t,s])) # Survivors
      }
    }
  }
  
  #-------------------------------#
  # 1.4) Harvest season (Oct-Mar) #
  #-------------------------------#
  
  for(s in 1:2){
    for(t in 1:Tmax){
      for(z in 1:Z){
        marN[z,t+1,s] ~ dbin(sH[z,t,s], octN_g[z,t,s]) # Survivors
        H[z,t+1,s] <- octN_g[z,t,s] - marN[z,t+1,s] # Harvests
      }
    }
  }

  #H[1:Z,1,1] <- 0
  H[1,1,1] <- 0
  H[2,1,1] <- 0
  H[3,1,1] <- 0
  #H[1:Z,1,2] <- 0
  H[1,1,2] <- 0
  H[2,1,2] <- 0
  H[3,1,2] <- 0
  
  #############################
  # 2) SIZE-AT-HARVEST MODULE #
  #############################
  
  # Data:
  # SaH[z,t] = Size-at-Harvest matrix (number of size class z individuals reported as harvested in the Oct[t-1] to Feb[t] harvest season)
  
  # Parameters:
  # H[z,t] = number of size class z individuals harvested in the Oct[t-1] to Feb[t] harvest season
  # ll[z,t] = recovery rate of size class z individuals harvested in the Oct[t-1] to Feb[t] harvest season
  
  for(s in 1:2){
    for(t in 2:Tmax){
      for(z in 1:Z){
        SaH[z,t,s] ~ dbin(ll[z,t,s], H[z,t,s])
      }
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
  
  for(s in 1:2){
    for(t in 1:Tmax){
      
      # Alive transitions (stochastic)
      
      ps[1,1,t,s] <- (1-gP[1,t,s])*sN[1,t,s]*sH[1,t,s] 
      ps[1,2,t,s] <- gP[1,t,s]*(1-gSL[t,s])*sN[2,t,s]*sH[2,t,s]
      ps[1,3,t,s] <- gP[1,t,s]*gSL[t,s]*sN[3,t,s]*sH[3,t,s]
      ps[1,4,t,s] <- (1-gP[1,t,s])*sN[1,t,s]*(1-sH[1,t,s])
      ps[1,5,t,s] <- gP[1,t,s]*(1-gSL[t,s])*sN[2,t,s]*(1-sH[2,t,s]) 
      ps[1,6,t,s] <- gP[1,t,s]*gSL[t,s]*sN[3,t,s]*(1-sH[3,t,s]) 
      ps[1,7,t,s] <- 1-sum(ps[1,1:6,t,s])
      
      ps[2,1,t,s] <- 0 
      ps[2,2,t,s] <- (1-gP[2,t,s])*sN[2,t,s]*sH[2,t,s]
      ps[2,3,t,s] <- gP[2,t,s]*sN[3,t,s]*sH[3,t,s] 
      ps[2,4,t,s] <- 0
      ps[2,5,t,s] <- (1-gP[2,t,s])*sN[2,t,s]*(1-sH[2,t,s]) 
      ps[2,6,t,s] <- gP[2,t,s]*sN[3,t,s]*(1-sH[3,t,s]) 
      ps[2,7,t,s] <- 1-sum(ps[2,1:6,t,s])
      
      ps[3,1,t,s] <- 0 
      ps[3,2,t,s] <- 0
      ps[3,3,t,s] <- sN[3,t,s]*sH[3,t,s]
      ps[3,4,t,s] <- 0 
      ps[3,5,t,s] <- 0
      ps[3,6,t,s] <- sN[3,t,s]*(1-sH[3,t,s])
      ps[3,7,t,s] <- 1-sum(ps[3,1:6,t,s])
      
      
      # Dead transitions (deterministic)
      
      #ps[4,1:6,t,s] <- 0
      ps[4,1,t,s] <- 0
      ps[4,2,t,s] <- 0
      ps[4,3,t,s] <- 0
      ps[4,4,t,s] <- 0
      ps[4,5,t,s] <- 0
      ps[4,6,t,s] <- 0
      ps[4,7,t,s] <- 1
      
      #ps[5,1:6,t,s] <- 0
      ps[5,1,t,s] <- 0
      ps[5,2,t,s] <- 0
      ps[5,3,t,s] <- 0
      ps[5,4,t,s] <- 0
      ps[5,5,t,s] <- 0
      ps[5,6,t,s] <- 0
      ps[5,7,t,s] <- 1
      
      #ps[6,1:6,t,s] <- 0
      ps[6,1,t,s] <- 0
      ps[6,2,t,s] <- 0
      ps[6,3,t,s] <- 0
      ps[6,4,t,s] <- 0
      ps[6,5,t,s] <- 0
      ps[6,6,t,s] <- 0
      ps[6,7,t,s] <- 1
      
      #ps[7,1:6,t,s] <- 0
      ps[7,1,t,s] <- 0
      ps[7,2,t,s] <- 0
      ps[7,3,t,s] <- 0
      ps[7,4,t,s] <- 0
      ps[7,5,t,s] <- 0
      ps[7,6,t,s] <- 0
      ps[7,7,t,s] <- 1
      
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
    
    for(t in 1:Tmax){
      
      po[1,1,t,s] <- pp[1,t,s]
      #po[1,2:6,t,s] <- 0
      po[1,2,t,s] <- 0
      po[1,3,t,s] <- 0
      po[1,4,t,s] <- 0
      po[1,5,t,s] <- 0
      po[1,6,t,s] <- 0
      po[1,7,t,s] <- 1-pp[1,t,s]
      
      po[2,1,t,s] <- 0
      po[2,2,t,s] <- pp[2,t,s]
      #po[2,3:6,t,s] <- 0
      po[2,3,t,s] <- 0
      po[2,4,t,s] <- 0
      po[2,5,t,s] <- 0
      po[2,6,t,s] <- 0
      po[2,7,t,s] <- 1-pp[2,t,s]
      
      #po[3,1:2,t,s] <- 0
      po[3,1,t,s] <- 0
      po[3,2,t,s] <- 0
      po[3,3,t,s] <- pp[3,t,s]
      #po[3,4:6,t,s] <- 0
      po[3,4,t,s] <- 0
      po[3,5,t,s] <- 0
      po[3,6,t,s] <- 0
      po[3,7,t,s] <- 1-pp[3,t,s]
      
      #po[4,1:3,t,s] <- 0
      po[4,1,t,s] <- 0
      po[4,2,t,s] <- 0
      po[4,3,t,s] <- 0
      po[4,4,t,s] <- ll[1,t,s]
      #po[4,5:6,t,s] <- 0
      po[4,5,t,s] <- 0
      po[4,6,t,s] <- 0
      po[4,7,t,s] <- 1-ll[1,t,s]
      
      #po[5,1:4,t,s] <- 0
      po[5,1,t,s] <- 0
      po[5,2,t,s] <- 0
      po[5,3,t,s] <- 0
      po[5,4,t,s] <- 0
      po[5,5,t,s] <- ll[2,t,s]
      po[5,6,t,s] <- 0
      po[5,7,t,s] <- 1-ll[2,t,s]
      
      #po[6,1:5,t,s] <- 0
      po[6,1,t,s] <- 0
      po[6,2,t,s] <- 0
      po[6,3,t,s] <- 0
      po[6,4,t,s] <- 0
      po[6,5,t,s] <- 0
      po[6,6,t,s] <- ll[3,t,s]
      po[6,7,t,s] <- 1-ll[3,t,s]
      
      #po[7,1:6,t,s] <- 0
      po[7,1,t,s] <- 0
      po[7,2,t,s] <- 0
      po[7,3,t,s] <- 0
      po[7,4,t,s] <- 0
      po[7,5,t,s] <- 0
      po[7,6,t,s] <- 0
      po[7,7,t,s] <- 1
      
    }
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
  		x[i, t] ~ dcat(ps[x[i, t-1], 1:7, t-1, CH.sex[i]])
  		
  		# Observation process: draw y(t) given x(t)
  		y[i, t] ~ dcat(po[x[i, t], 1:7, t, CH.sex[i]])
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
  
  for(s in 1:2){
    for(z in 1:Z){
      
      sN[z,1:Tmax,s] <- exp(-mN[z,1:Tmax,s])
      sH[z,1:Tmax,s] <- exp(-mH[z,1:Tmax,s])
      
      log(mN[z,1:Tmax,s]) <- log(Mu.mN[z,s]) + epsilon.mN[1:Tmax,s]
      log(mH[z,1:Tmax,s]) <- log(Mu.mH[z,s]) + epsilon.mH[1:Tmax,s]
      
      Mu.mN[z,s] ~ dunif(0, 5)
      Mu.mH[z,s] ~ dunif(0, 5)
      
    }
    
    for(t in 1:Tmax){
      
      epsilon.mN[t,s] ~ dnorm(0, sd = sigma.mN[s])
      epsilon.mH[t,s] ~ dnorm(0, sd = sigma.mH[s])
    }
    
    sigma.mN[s] ~ dunif(0, 5)
    sigma.mH[s] ~ dunif(0, 5)
  }
  
  
  #------------------------#
  # 5.2) Growth parameters #
  #------------------------#
  
  for(s in 1:2){
    for(z in 1:(Z-1)){
      
      logit(gP[z,1:Tmax,s]) <- logit(Mu.gP[z,s]) + epsilon.gP[z,1:Tmax,s]
      
      Mu.gP[z,s] ~ dunif(0, 1)
      sigma.gP[z,s] ~ dunif(0, 5)
    }
    
    
    for(t in 1:Tmax){
      
      for(z in 1:(Z-1)){
        epsilon.gP[z,t,s] ~ dnorm(0, sd = sigma.gP[z,s])
      }
      
      epsilon.gSL[t,s] ~ dnorm(0, sd = sigma.gSL[s])
    }
    
    logit(gSL[1:Tmax,s]) <- logit(Mu.gSL[s]) + epsilon.gSL[1:Tmax,s]
    
    Mu.gSL[s] ~ dunif(0, 1)
    sigma.gSL[s] ~ dunif(0, 5)
  }
  
  
  #------------------------------------#
  # 5.3) Recapture/Recovery parameters #
  #------------------------------------#
  
  for(s in 1:2){
    for(z in 1:Z){
      
      logit(pp[z,1:Tmax,s]) <- logit(Mu.pp[z,s]) + epsilon.pp[1:Tmax,s]
      
      logit(ll[z,1:Tmax,s]) <- logit(Mu.ll[s]) + epsilon.ll[1:Tmax,s]
      
      Mu.pp[z,s] ~ dunif(0, 1)
    }
    
    Mu.ll[s] ~ dunif(0, 1)
    
    
    for(t in 1:Tmax){
      
      epsilon.pp[t,s] ~ dnorm(0, sd = sigma.pp[s])
      epsilon.ll[t,s] ~ dnorm(0, sd = sigma.ll[s])
    }
    
    sigma.pp[s] ~ dunif(0, 5)
    sigma.ll[s] ~ dunif(0, 5)
  }

  
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
  for(s in 1:2){
    for(z in 1:Z){
      initN[z,s] ~ dlnorm(3.5, 1)
      marN[z,1,s] <- round(initN[z,s])
    }
  }

  
  
  #----------------------------------------#
  # 5.5) Early survival ('Free' Parameter) #
  #----------------------------------------#
  
  for(s in 1:2){
    for(t in 1:Tmax){
      
      s0[t,s] <- exp(-m0[t,s])
      log(m0[t,s]) <- log(Mu.m0[s]) + epsilon.m0[t,s]
      epsilon.m0[t,s] ~ dnorm(0, sd = sigma.m0[s])
    }
    
    Mu.m0[s] ~ dunif(0, 5)
    sigma.m0[s] ~ dunif(0, 5)
  }
  
  
})


#*******************************************************************************#
#* INITIAL VALUES
#*******************************************************************************#

## Load functions for simulating initial values
source('WildBoarIPM_InitialValuesSim_TwoSex.R')


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
				#'Off', 
				'Off_tot', 'YOY'
				)

## MCMC settings
ni <- 100000
nb <- 50000
nt <- 10
nc <- 2

## Sample initial values (2 chains)
Inits <- list(
  inits.duplicate(WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5))), 
  inits.duplicate(WB.IPM.inits.convert(model = 'A', const.N1 = c(rowSums(SaH)), const.marProps = diag(3), extra.N = round(c(rowSums(SaH))/5))))


## Reformat initial values for unlinked model
# for(i in 1:nc){
#   Inits[[i]]$Off_tot <- abind::abind(Inits[[i]]$Off_tot/2, Inits[[i]]$Off_tot/2, along = 3)
#   Inits[[i]]$Off <- NULL
# }

#*******************************************************************************#
#* RUN
#*******************************************************************************#

t1 <- Sys.time()
WB.IPM.TwoSex <- nimbleMCMC(code = WB.code, constants = WB.constants, data = WB.data, inits = Inits, monitors = parameters, niter = ni, nburnin = nb, nchains = nc, thin = nt, setSeed = mySeed, samplesAsCodaMCMC = TRUE)
t2 <- Sys.time()
t2-t1

saveRDS(WB.IPM.TwoSex, file = 'WildBoarIPM_TwoSex_MCMCsamples.rds')

pdf('WildBoarIPM_TwoSex_traces.pdf', height = 8, width = 11)
plot(WB.IPM.TwoSex)
dev.off()

WB.IPM.TwoSex.trim <- mcmc.list(
  as.mcmc(WB.IPM.TwoSex[[1]][2500:5000, ]),
  as.mcmc(WB.IPM.TwoSex[[2]][2500:5000, ])
)

#*******************************************************************************#
#* COMPARE FEMALE VS. MALE PARAMETERS 
#*******************************************************************************#

library(reshape2)
library(ggplot2)

## Set up parameter matching across sexes 
param.match <- data.frame(
  Parameter = c(paste0('Mu.mN[',1:3,', 1]'), paste0('Mu.mH[',1:3,', 1]'), 'sigma.mN[1]', 'sigma.mH[1]', 
                paste0('Mu.gP[',1:2,', 1]'), 'Mu.gSL[1]', 'sigma.gP[1]', 'sigma.gSL[1]', 
                paste0('Mu.pp[',1:3,', 1]'), 'Mu.ll[1]', 'sigma.pp[1]', 'sigma.ll[1]',
                'Mu.m0[1]', 'sigma.m0[1]',
                paste0('marN[1, ', 1:27, ', 1]'), paste0('marN[2, ', 1:27, ', 1]'), paste0('marN[3, ', 1:27, ', 1]'),
                paste0('H[1, ', 1:26, ', 1]'), paste0('H[2, ', 1:26, ', 1]'), paste0('H[3, ', 1:26, ']'),
                paste0('Off[', 1:26, ', 1]'), 
                paste0('YOY[', 1:26, ', 1]'), 
                
                paste0('Mu.mN[',1:3,', 2]'), paste0('Mu.mH[',1:3,', 2]'), 'sigma.mN[2]', 'sigma.mH[2]', 
                paste0('Mu.gP[',1:2,', 2]'), 'Mu.gSL[2]', 'sigma.gP[2]', 'sigma.gSL[2]', 
                paste0('Mu.pp[',1:3,', 2]'), 'Mu.ll[2]', 'sigma.pp[2]', 'sigma.ll[2]',
                'Mu.m0[2]', 'sigma.m0[2]',
                paste0('marN[1, ', 1:27, ', 2]'), paste0('marN[2, ', 1:27, ', 2]'), paste0('marN[3, ', 1:27, ', 2]'),
                paste0('H[1, ', 1:26, ', 2]'), paste0('H[2, ', 1:26, ', 2]'), paste0('H[3, ', 1:26, ', 2]'),
                paste0('Off[', 1:26, ', 2]'), 
                paste0('YOY[', 1:26, ', 2]')
  ),
  ParamGeneral = rep(c(paste0('Mu.mN[',1:3,']'), paste0('Mu.mH[',1:3,']'), 'sigma.mN', 'sigma.mH', 
                       paste0('Mu.gP[',1:2,']'), 'Mu.gSL', 'sigma.gP', 'sigma.gSL', 
                       paste0('Mu.pp[',1:3,']'), 'Mu.ll', 'sigma.pp', 'sigma.ll',
                       'Mu.m0', 'sigma.m0',
                       paste0('marN[1, ', 1:27, ']'), paste0('marN[2, ', 1:27, ']'), paste0('marN[3, ', 1:27, ']'),
                       paste0('H[1, ', 1:26, ']'), paste0('H[2, ', 1:26, ']'), paste0('H[3, ', 1:26, ']'),
                       paste0('Off[', 1:26, ']'),  
                       paste0('YOY[', 1:26, ']')), 2),
  Sex = rep(c('F', 'M'), each = 232)
)

## Make subsets of posterior data for different parameters
post.data <- melt(as.matrix(WB.IPM.TwoSex.trim))
colnames(post.data) <- c('Sample', 'Parameter', 'Value')

post.data <- merge(post.data, param.match, by = 'Parameter', all.x = T)

sub.data <- subset(post.data, !is.na(Sex) & !is.na(ParamGeneral))

VR.data <- subset(sub.data, ParamGeneral %in% c(paste0('Mu.mN[',1:3,']'), paste0('Mu.mH[',1:3,']'), 'sigma.mN', 'sigma.mH', 
                                                paste0('Mu.gP[',1:2,']'), 'Mu.gSL', 'sigma.gP', 'sigma.gSL', 
                                                paste0('Mu.pp[',1:3,']'), 'Mu.ll', 'sigma.pp', 'sigma.ll',
                                                'Mu.m0', 'sigma.m0'))
N.data <- subset(sub.data, ParamGeneral %in% c(paste0('marN[1, ', 1:27, ']'), paste0('marN[2, ', 1:27, ']'), paste0('marN[3, ', 1:27, ']')))
N.data$ParamGeneral <- factor(N.data$ParamGeneral, levels = c(paste0('marN[1, ', 1:27, ']'), paste0('marN[2, ', 1:27, ']'), paste0('marN[3, ', 1:27, ']')))

H.data <- subset(sub.data, ParamGeneral %in% c(paste0('H[1, ', 1:27, ']'), paste0('H[2, ', 1:27, ']'), paste0('H[3, ', 1:27, ']')))
H.data$ParamGeneral <- factor(H.data$ParamGeneral, levels = c(paste0('H[1, ', 1:27, ']'), paste0('H[2, ', 1:27, ']'), paste0('H[3, ', 1:27, ']')))

Off.data <- subset(sub.data, ParamGeneral %in% c(paste0('Off[', 1:26, ']')))
Off.data$ParamGeneral <- factor(Off.data$ParamGeneral, levels = c(paste0('Off[', 1:26, ']')))

YOY.data <- subset(sub.data, ParamGeneral %in% c(paste0('YOY[', 1:26, ']')))
YOY.data$ParamGeneral <- factor(YOY.data$ParamGeneral, levels = c(paste0('YOY[', 1:26, ']')))

## Plot overlap between male and female parameters

# Vital rates
ggplot(VR.data, aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Population size (class 1)
ggplot(subset(N.data, ParamGeneral%in%c(paste0('marN[1, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Population size (class 2)
ggplot(subset(N.data, ParamGeneral%in%c(paste0('marN[2, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Population size (class 3)
ggplot(subset(N.data, ParamGeneral%in%c(paste0('marN[3, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Harvest (class 1)
ggplot(subset(H.data, ParamGeneral%in%c(paste0('H[1, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Harvest (class 2)
ggplot(subset(H.data, ParamGeneral%in%c(paste0('H[2, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Harvest (class 3)
ggplot(subset(H.data, ParamGeneral%in%c(paste0('H[3, ', 1:27, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# Offspring 
ggplot(subset(Off.data, ParamGeneral%in%c(paste0('Off[', 1:26, ']'))), aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()

# YOY
ggplot(YOY.data, aes(x = Value, group = Sex)) +
  geom_density(aes(color = Sex, fill = Sex), alpha = 0.5) +
  facet_wrap(~ ParamGeneral, scale = 'free') +
  theme_bw()
