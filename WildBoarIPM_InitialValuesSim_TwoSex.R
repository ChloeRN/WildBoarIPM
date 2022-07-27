#*******************************************************************************#
#* INITIAL VALUE FUNCTIONS
#*******************************************************************************#

## Function to set initial values for latent states (mark-recapture-recovery model)
state.inits <- function(ch, f){
  
  ## Turn unknown states into NA
  ch[ch==7] <- NA
  
  ## States after newly dead from harvest become "dead"
  ## --> Individuals in states 4, 5, and 6 always transition to state 7 
  vH <- which(ch==4 | ch==5 | ch==6, arr.ind = T)
  for (i in 1:nrow(vH)){
    ifelse(vH[i,2]!=ncol(ch), ch[vH[i,1], (vH[i,2]+1):ncol(ch)] <- 7, next)}
  
  ## Remaining states
  for(i in 1:nrow(ch)){
    for(x in (f[i]+1):dim(ch)[2]){
      if(is.na(ch[i,x])){
        
        ObsList <- which(!is.na(ch[i,x:dim(ch)[2]]))
        nextObsIndex <- ifelse(length(ObsList) > 0, min(ObsList) + x - 1, NA)
        nextObs <- ifelse(!is.na(nextObsIndex), ch[i, nextObsIndex], NA)
        
        if(ch[i,x-1]==1){ # If the previous capture was in state 1...
          
          if(!is.na(nextObs)){
            if(nextObs%in%c(1,4)){ch[i,x] <- 1}else{ch[i,x] <- 2}
          }else{ch[i,x] <- 2}
          
          #... then next unknown LS is 1 if next capture is in state 1 or 4, and 2 otherwise        	
        }
        
        if(ch[i,x-1]==2){ # If the previous capture was in state 2...
          
          if(!is.na(nextObs)){
            if(nextObs%in%c(2,5)){ch[i,x] <- 2}else{ch[i,x] <- 3}
          }else{ch[i,x] <- 3}
          
          #... then next unknown LS is 2 if next capture is in state 2 or 5, and 3 otherwise        	
        }
        
        if(ch[i,x-1]==3){ # If the previous capture was in state 3...
          ch[i,x] <- 3
          #... then next unknown LS is 3       	
        }        
      }
    }
  }
  return(ch)
}


## Function to simulate and assemble intial values (multinomial model version)

#  NOTE: The function 'WB.IPM.inits()' was built for a generalized version of
#        the model that models growth (state transitions) as multinomial random
#        variables. The following function 'WB.IPM.inits.convert()' converts the
#        initial values to match to the model we fit (sequential binomial formulation)

WB.IPM.inits <- function(model, const.N1, const.marProps, extra.N) {
  Z <- WB.constants$Z
  Tmax <- WB.constants$Tmax+1
  
  # Initial values
  Mu.ll <- 1
  sigma.ll <- runif(1, 0, 2)
  epsilon.ll <- rep(0, WB.constants$Tmax)
  
  Mu.pp <- runif(3, 0.2, 0.7)
  sigma.pp <- runif(1, 0, 2)
  epsilon.pp <- rep(0, WB.constants$Tmax)
  
  Mu.mH <- c(runif(1, 0.2, 0.5), runif(1, 0.2, 0.5), runif(1, 0.2, 0.5))
  sigma.mH <- runif(1, 0, 2)
  epsilon.mH <- rep(0, WB.constants$Tmax)
  #epsilon.mH <- rnorm(WB.constants$Tmax, 0, sigma.mH)
  
  Mu.pB <- matrix(runif(9, 0.8, 1), nrow = 3)
  Mu.nF <- runif(3, 8, 10)
  sigma.pB <- runif(1, 0, 2)
  sigma.nF <- runif(1, 0, 2)
  epsilon.pB <- rep(0, WB.constants$Tmax)
  epsilon.nF <- rep(0, WB.constants$Tmax)
  #epsilon.pB <- rnorm(WB.constants$Tmax, 0, sigma.pB)
  #epsilon.nF <- rnorm(WB.constants$Tmax, 0, sigma.nF)
  
  Mu.m0 <- runif(1, 0.1, 0.2)
  sigma.m0 <- runif(1, 0, 1)
  epsilon.m0 <- rep(0, WB.constants$Tmax)
  #epsilon.m0 <- rnorm(WB.constants$Tmax, 0, sigma.m0)
  
  Mu.mN <- runif(3, 0.1, 0.5)
  sigma.mN <- runif(1, 0, 2)
  epsilon.mN <- rep(0, WB.constants$Tmax)
  #epsilon.mN <- rnorm(WB.constants$Tmax, 0, sigma.mN)
  
  G <- array(0, dim = c(2, 3, Tmax))
  for (tt in 1:Tmax) {
    
    #G[1, 1:3, tt] <- rdirichlet(1, c(1, 1, 1))
    #G[2, 2:3, tt] <- rdirichlet(1, c(1, 1))
    
    G[1, 1, tt] <- runif(1, 0, 0.01)
    G[1, 2, tt] <- runif(1, 0.1, 0.9)
    G[1, 3, tt] <- 1 - G[1, 2, tt] - G[1, 1, tt]
    
    G[2, 3, tt] <- runif(1, 0.6, 0.9)
    G[2, 2, tt] <- 1 - G[2, 3, tt]
  }
  
  # Empty vectors and matrices
  H <- matrix(0, nrow = Z, ncol = Tmax)
  marN <- array(0, dim = c(Z, 2 * Z, Tmax))
  octN <- array(0, dim = c(Z, 2 * Z, Tmax))
  N <- matrix(NA, nrow = Z, ncol = Tmax)
  Off <- matrix(NA, nrow = Z, ncol = Tmax-1)
  YOY <- YOY.all <- rep(NA, Tmax-1)
  
  N[1:Z, 1] <- c(const.N1[1], const.N1[2], const.N1[3])
  
  
  for (t in 2:Tmax) {
    
    # a) Set number of individuals harvested
    #---------------------------------------
    
    ## a.1) Assume true number of individuals harvested = number reported
    for (z in 1:Z) {
      if(t < Tmax){
        H[z, t] <- SaH[z, t, 1]
      }else{
        H[z, t] <- round(sum(SaH[z,,1])/dim(SaH)[2])
        # --> estimating as mean since no data available
      }
      
    }
    
    # b) Set March census population sizes 
    #-------------------------------------
    
    ## b.1) Set summarised census population sizes (N) to the specified constant values
    N[1:Z, t] <- c(const.N1[1], const.N1[2], const.N1[3])
    
    ## b.2) Fill harvests into detailed census population sizes (marN)
    for(z in 1:z){
      marN[z, z+3, t] <- H[z, t] # Note: this only works for models A and B
    }
    
    ## b.3) Fill survivors into detailed census population sizes (marN), split by from-class
    for(z in 1:Z){
      splitN <- rmultinom(1, size = N[z, t], prob = const.marProps[z,])
      marN[1:Z, z, t] <- splitN
    }
    
    # c) Set October population sizes
    #--------------------------------
    
    ## c.1) Calculate source state (= from-class) proportions for survivors (Mar - Sep)
    octPropsA <- matrix(NA, nrow = Z, ncol = Z)
    
    octPropsA[1, ] <- c(1, 0, 0)
    # --> State 1 survivors can only have been produced by state 1 individuals (because there is no shrinkage)
    
    mN <- exp(log(Mu.mN[1:Z]) + epsilon.mN[t])
    sN <- exp(-mN[1:Z])
    
    if (model == 'A') {
      octPropsA[2, ] <-
        c(G[1, 2, t], G[2, 2, t], 0) / (G[1, 2, t] + G[2, 2, t])
      octPropsA[3, ] <-
        c(G[1, 3, t], G[2, 3, t], 1) / (G[1, 3, t] + G[2, 3, t] + 1)
    }
    
    if (model == 'B') {
      octPropsA[2, ] <-
        c(0, sN[1] * G[1, 2, t], sN[2] * G[2, 2, t]) / (sN[1] * G[1, 2, t] + sN[2] *
                                                          G[2, 2, t])
      octPropsA[3, ] <-
        c(sN[1] * G[1, 3, t], sN[2] * G[2, 3, t], sN[3]) / (sN[1] * G[1, 3, t] + sN[2] *
                                                              G[2, 3, t] + sN[3])
    }
    
    ## c.2) Fill survivors into detailed census population sizes (marN), split by from-class
    for(z in 1:Z){
      octN[1:Z, z, t-1] <- rmultinom(1, size = sum(marN[z,,t]), prob = octPropsA[z,])
    }
    
    ## c.3) Calculate source state (= from-class) proportions of natural deaths (Mar - Sep)
    octPropsD <- matrix(NA, nrow = Z, ncol = Z)
    
    octPropsD[1, ] <- c(1, 0, 0)
    
    
    if (model == 'A') {
      octPropsD[2, ] <-
        c(G[1, 2, t], G[2, 2, t], 0) / (G[1, 2, t] + G[2, 2, t])
      octPropsD[3, ] <-
        c(G[1, 3, t], G[2, 3, t], 1) / (G[1, 3, t] + G[2, 3, t] + 1)
    }
    
    if (model == 'B') {
      octPropsD[2, ] <-
        c(0, (1 - sN[1]) * G[1, 2, t], (1 - sN[2]) * G[2, 2, t]) / ((1 - sN[1]) *
                                                                      G[1, 2, t] + (1 - sN[2]) * G[2, 2, t])
      octPropsD[3, ] <-
        c((1 - sN[1]) * G[1, 3, t], (1 - sN[2]) * G[2, 3, t], (1 - sN[3])) / ((1 -
                                                                                 sN[1]) * G[1, 3, t] + (1 - sN[2]) * G[2, 3, t] + (1 - sN[3]))
    }
    
    ## c.4) Add a number of "extra" individuals that died (Mar - Sep)
    for(z in 1:Z){
      octN[1:Z, z+3, t-1] <- rmultinom(1, size = extra.N[z], prob = octPropsD[z,])
    }
    
    # d) Quantify reproduction and surplus individuals
    #-------------------------------------------------
    
    ## d.1) Count the number of discrepant (missing/surplus) individual in each class
    discrep.N <- rowSums(octN[,,t-1]) - N[,t]
    
    ## d.2) Add missing individuals as natural deaths
    for(z in 1:Z){
      if(discrep.N[z] < 0){
        octN[z,z+3,t-1] <- octN[z,z+3,t-1] + abs(discrep.N[z])	
      }
    }
    # NOTE: I am adding them to the non-growing fraction of each class (but they could also be distributed)
    
    ## d.3) Remove surplus deaths (for from-classes 2 and 3)
    for(z in 2:Z){
      if(discrep.N[z] > 0){
        octN[z,z+3,t-1] <- octN[z,z+3,t-1] - abs(discrep.N[z])	
      }
    }
    
    ## d.4) Extract the total number of YOY (class 1 surplus)
    YOY[t-1] <- discrep.N[1]
    
    ## d.5) Add offspring that died before recruiting
    YOY.all[t-1] <- round(YOY[t-1]*(2-exp(-Mu.m0)))
    # NOTE: Not really necessary
    
    ## d.6) Distribute YOY across mother classes
    matProps <- (N[,t-1]*Mu.pB[,2]*Mu.nF)/sum(N[,t-1]*Mu.pB[,2]*Mu.nF)
    for(z in 1:Z){
      Off[1:Z, t-1] <- rmultinom(1, size = YOY.all[t-1], prob = matProps[1:Z])
    }
    
  }
  
  return(list(
    N = N, 
    initN = N[,1],
    octN = octN, 
    marN = marN, 
    H = H, 
    Off = Off, 
    YOY = YOY, 
    Mu.ll = Mu.ll, sigma.ll = sigma.ll, epsilon.ll = epsilon.ll,
    Mu.pp = Mu.pp, sigma.pp = sigma.pp, epsilon.pp = epsilon.pp,
    Mu.mH = Mu.mH, sigma.mH = sigma.mH, epsilon.mH = epsilon.mH,
    Mu.pB = Mu.pB, sigma.pB = sigma.pB, epsilon.pB = epsilon.pB, 
    Mu.nF = Mu.nF, sigma.nF = sigma.nF, epsilon.nF = epsilon.nF,
    Mu.m0 = Mu.m0, sigma.m0 = sigma.m0, epsilon.m0 = epsilon.m0,
    Mu.mN = Mu.mN, sigma.mN = sigma.mN, epsilon.mN = epsilon.mN,	
    G = G,
    x = state.inits(mydata, first)
  ))
  
}


## Function to convert initial values for multinomial model to initial values for sequential multinomial model

WB.IPM.inits.convert <- function(model, const.N1, const.marProps, extra.N) {
  
  Z <- WB.constants$Z
  Tmax <- WB.constants$Tmax+1
  
  ## Simulate initial values for multinomial model
  origI <- WB.IPM.inits(model, const.N1, const.marProps, extra.N)
  
  ## Make empty matrices/vectors
  marN_g <- array(NA, dim = c(Z, Z, Tmax-1))
  marN1_plus <- rep(NA, Tmax)
  octN_g <- matrix(NA, nrow = Z, ncol = Tmax-1)
  gP <- matrix(NA, nrow = Z-1, ncol = Tmax)
  gSL <- rep(NA, Tmax)
  Mu.gP <- rep(NA, Z-1)
  
  ## Convert March population census values
  marN <- origI$N
  
  for(t in 1:(Tmax-1)){
    
    ## Extract post-growth march census values
    marN_g[1:Z, 1:Z, t] <- origI$octN[1:Z, 1:Z, t] + origI$octN[1:Z, (1:Z)+Z, t]
    marN1_plus[t] <- sum(marN_g[1, 2:3, t])
    
    ## Summarize October census numbers
    for(z in 1:Z){		
      octN_g[z, t] <- sum(origI$octN[1:Z, z, t])
    }
  }
  
  for(t in 1:Tmax){
    
    ## Extract time-dependent binomial growth parameters
    gP[1, t] <- sum(origI$G[1, 2:3, t])
    gP[2, t] <- origI$G[2, 3, t]
    gSL[t] <- origI$G[1, 3, t]/gP[1, t]
  }
  
  ## Calculate time averages for binomial growth parameters
  for(z in 1:(Z-1)){
    Mu.gP[z] <- mean(gP[z,])
  }
  Mu.gSL <- mean(gSL)
  
  ## Sample values for additional parameters
  sigma.gP <- runif(2, 0, 2)
  sigma.gSL <- runif(1, 0, 2)
  epsilon.gP <- matrix(0,nrow = Z-1, ncol = Tmax)
  epsilon.gSL <- rep(0, Tmax)
  
  ## Assemble converted initial values
  return(list(
    marN = marN, 
    initN = marN[,1],
    marN_g = marN_g,
    marN1_plus = marN1_plus,
    octN_g = octN_g,
    H = origI$H, 
    Off = origI$Off, 
    YOY = origI$YOY, 
    Mu.ll = origI$Mu.ll, sigma.ll = origI$sigma.ll, epsilon.ll = origI$epsilon.ll,
    Mu.pp = origI$Mu.pp, sigma.pp = origI$sigma.pp, epsilon.pp = origI$epsilon.pp,
    Mu.mH = origI$Mu.mH, sigma.mH = origI$sigma.mH, epsilon.mH = origI$epsilon.mH,
    Mu.pB = origI$Mu.pB, sigma.pB = origI$sigma.pB, epsilon.pB = origI$epsilon.pB, 
    Mu.nF = origI$Mu.nF, sigma.nF = origI$sigma.nF, epsilon.nF = origI$epsilon.nF,
    Mu.m0 = origI$Mu.m0, sigma.m0 = origI$sigma.m0, epsilon.m0 = origI$epsilon.m0,
    Mu.mN = origI$Mu.mN, sigma.mN = origI$sigma.mN, epsilon.mN = origI$epsilon.mN,	
    Mu.gP = Mu.gP, sigma.gP = sigma.gP, epsilon.gP = epsilon.gP,
    Mu.gSL = Mu.gSL, sigma.gSL = sigma.gSL, epsilon.gSL = epsilon.gSL,
    x = state.inits(mydata, first)
  ))
}

## Function to duplicate initial values for assigning identical inital values to 
## female and male parameters

## NOTE: This only works for the test situation which duplicates female data to 
##       "simulate" male data. 
##       For including actual male data, the initial value simulation functions
##       above will need to be updated.

inits.duplicate <- function(Inits){
  
  Off_tot <- Inits$Off*2
  Off <- cbind(colSums(Inits$Off), colSums(Inits$Off))
  
  return(list(
    marN = abind::abind(Inits$marN, Inits$marN, along = 3), 
    initN = cbind(Inits$marN[,1], Inits$marN[,1]),
    marN_g = abind::abind(Inits$marN_g, Inits$marN_g, along = 4),
    marN1_plus = cbind(Inits$marN1_plus, Inits$marN1_plus),
    octN_g = abind::abind(Inits$octN_g, Inits$octN_g, along = 3),
    H = abind::abind(Inits$H, Inits$H, along = 3), 
    Off_tot = Off_tot, 
    Off = Off,
    YOY = cbind(Inits$YOY, Inits$YOY), 
    Mu.ll = rep(Inits$Mu.ll, 2), sigma.ll = rep(Inits$sigma.ll, 2), epsilon.ll = cbind(Inits$epsilon.ll, Inits$epsilon.ll),
    Mu.pp = cbind(Inits$Mu.pp, Inits$Mu.pp), sigma.pp = rep(Inits$sigma.pp, 2), epsilon.pp = cbind(Inits$epsilon.pp, Inits$epsilon.pp),
    Mu.mH = cbind(Inits$Mu.mH, Inits$Mu.mH), sigma.mH = rep(Inits$sigma.mH, 2), epsilon.mH = cbind(Inits$epsilon.mH, Inits$epsilon.mH),
    Mu.pB = Inits$Mu.pB, sigma.pB = Inits$sigma.pB, epsilon.pB = Inits$epsilon.pB, 
    Mu.nF = Inits$Mu.nF, sigma.nF = Inits$sigma.nF, epsilon.nF = Inits$epsilon.nF,
    Mu.m0 = rep(Inits$Mu.m0, 2), sigma.m0 = rep(Inits$sigma.m0, 2), epsilon.m0 = cbind(Inits$epsilon.m0, Inits$epsilon.m0),
    Mu.mN = cbind(Inits$Mu.mN, Inits$Mu.mN), sigma.mN = rep(Inits$sigma.mN, 2), epsilon.mN = cbind(Inits$epsilon.mN, Inits$epsilon.mN),	
    Mu.gP = cbind(Inits$Mu.gP, Inits$Mu.gP), sigma.gP = cbind(Inits$sigma.gP, Inits$Mu.gP), epsilon.gP = abind::abind(Inits$epsilon.gP,Inits$epsilon.gP, along = 3),
    Mu.gSL = rep(Inits$Mu.gSL, 2), sigma.gSL = rep(Inits$sigma.gSL, 2), epsilon.gSL = cbind(Inits$epsilon.gSL,Inits$epsilon.gSL),
    x = Inits$x
  ))
}
