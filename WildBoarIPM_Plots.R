## Load packages
library(mcmcplots)
library(coda)
library(data.table)
library(ggplot2)
library(viridis)

## Load IPM output
load('WildBoarIPM_MCMCsamples.RData')

out.mat <- as.matrix(WB.IPM)

###########################
#### EXTRACT ESTIMATES ####
###########################

## Set the relevant years
minY <- 1991
maxY <- 2016
years <- minY:maxY

## Prepare arrays for storing population size and vital rate estimates
marN <- octN <- H <- matrix(NA, nrow = length(years), ncol = 3, dimnames = list(years, c('lCI', 'Median', 'uCI')))
dN <- dH <- nF <- pB <- array(NA, dim = c(3, length(years), 3), dimnames = list(c('lCI', 'Median', 'uCI'), years, c('Small', 'Medium', 'Large')))
s0 <- gSM <- gSL <- gML <- matrix(NA, nrow = length(years), ncol = 3, dimnames = list(years, c('lCI', 'Median', 'uCI')))
marProp <- octProp <- HProp <- matrix(NA, nrow = length(years), ncol = 3, dimnames = list(years, c('Small', 'Medium', 'Large')))

## Extract posterior summaries for year-specific quantities
for(t in 1:length(years)){
  
  # Collect posterior samples
  sam.marN <- out.mat[,paste0('marN[1, ',t,']')] + out.mat[,paste0('marN[2, ',t,']')] + out.mat[,paste0('marN[3, ',t,']')]
  sam.octN <- out.mat[,paste0('octN_g[1, ',t,']')] + out.mat[,paste0('octN_g[2, ',t,']')] + out.mat[,paste0('octN_g[3, ',t,']')]
  sam.H <- out.mat[,paste0('H[1, ',t,']')] + out.mat[,paste0('H[2, ',t,']')] + out.mat[,paste0('H[3, ',t,']')]
  
  sam.s0 <- exp(-out.mat[,paste0('m0[',t,']')])
  sam.gSM <- out.mat[,paste0('gP[1, ',t,']')]*(1-out.mat[,paste0('gSL[',t,']')])
  sam.gSL <- out.mat[,paste0('gP[1, ',t,']')]*out.mat[,paste0('gSL[',t,']')]
  sam.gML <- out.mat[,paste0('gP[2, ',t,']')]

  for(z in 1:3){
    sam.dN <- 1 - exp(-out.mat[,paste0('mN[',z,', ',t,']')])
    sam.dH <- 1 - exp(-out.mat[,paste0('mH[',z,', ',t,']')])
    sam.nF <- out.mat[,paste0('nF[',z,', ',t,']')]
    sam.pB <- out.mat[,paste0('pB[',z,', ',t,']')]
    
    sam.marProp <- out.mat[,paste0('marN[',z,', ',t,']')]/sam.marN
    sam.octProp <- out.mat[,paste0('octN_g[',z,', ',t,']')]/sam.octN
    
    if(t > 1){sam.HProp <- out.mat[,paste0('H[',z,', ',t,']')]/sam.H}
    
    # Summarize samples
    dN[,t,z] <- quantile(sam.dN, probs = c(0.025, 0.5, 0.995))
    dH[,t,z] <- quantile(sam.dH, probs = c(0.025, 0.5, 0.995))
    nF[,t,z] <- quantile(sam.nF, probs = c(0.025, 0.5, 0.995))
    pB[,t,z] <- quantile(sam.pB, probs = c(0.025, 0.5, 0.995))
    
    marProp[t,z] <- quantile(sam.marProp, probs = 0.5)
    octProp[t,z] <- quantile(sam.octProp, probs = 0.5)
    if(t > 1){HProp[t,z] <- quantile(sam.HProp, probs = 0.5)}
  }

  marN[t,] <- quantile(sam.marN, probs = c(0.025, 0.5, 0.995))
  octN[t,] <- quantile(sam.octN, probs = c(0.025, 0.5, 0.995))
  if(t > 1){H[t,] <- quantile(sam.H, probs = c(0.025, 0.5, 0.995))}
  
  s0[t,] <- quantile(sam.s0, probs = c(0.025, 0.5, 0.995)) 
  gSM[t,] <- quantile(sam.gSM, probs = c(0.025, 0.5, 0.995)) 
  gSL[t,] <- quantile(sam.gSL, probs = c(0.025, 0.5, 0.995)) 
  gML[t,] <- quantile(sam.gML, probs = c(0.025, 0.5, 0.995)) 
}

##################
#### PLOTTING ####
##################

# 1) Total population sizes over time #
#-------------------------------------#

## Arrange data
marN.data <- as.data.frame(marN)
marN.data$Year <- years

octN.data <- as.data.frame(octN)
octN.data$Year <- years + 0.5

H.data <- as.data.frame(H)
H.data$Year <- years-0.25

C.S <- c(166, 94, 49, 85, 124, 198, 139, 116, 112, 177, 143, 170, 224, 173, 184, 186, 179, 265, 196, 164, 123, 177, 181, 119,  206, 148)
C.M <- c(56, 71, 31, 31, 75, 127, 53, 108, 100, 58, 68, 125, 99, 66, 78, 68, 187, 118, 64, 83, 108, 91, 31, 102, 69, 45)
C.L <- c(34, 34, 15, 20, 57, 90, 63, 51, 65, 39, 51, 46, 75, 31, 19, 36, 42, 38, 59, 42, 33, 32, 34, 31, 70, 33)
C.data <- data.frame(CountH = C.S + C.M + C.L, Year = years-0.25)

## Plot population size estimates
ggplot(marN.data) + geom_point(aes(x = Year, y = Median)) + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI)) + 
  geom_point(data = octN.data, aes(x = Year, y = Median), color = '#44A694') + 
  geom_linerange(data = octN.data, aes(x = Year, ymin = lCI, ymax = uCI), color = '#44A694') +
  geom_point(data = H.data, aes(x = Year, y = Median), color = '#7F4AA9') + 
  geom_linerange(data = H.data, aes(x = Year, ymin = lCI, ymax = uCI), color = '#7F4AA9') +
  geom_point(data = C.data, aes(x = Year, y = CountH), color = 'hotpink', shape = 3) +
  ylab('Number of females') + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

# 2) Natural and hunting mortality over time #
#--------------------------------------------#

## Arrange data
dN.data <- as.data.frame(rbind(t(dN[,,1]), t(dN[,,2]), t(dN[,,3])))
dN.data$Year <- rep(years, 3)
dN.data$Year2 <- c(years-0.2, years, years+0.2)
dN.data$SizeClass <- rep(c('Small', 'Medium', 'Large'), each = length(years))
rownames(dN.data) <- NULL

dH.data <- as.data.frame(rbind(t(dH[,,1]), t(dH[,,2]), t(dH[,,3])))
dH.data$Year <- rep(years, 3)
dH.data$Year2 <- c(years-0.2, years, years+0.2)
dH.data$SizeClass <- rep(c('Small', 'Medium', 'Large'), each = length(years))
rownames(dH.data) <- NULL

## Plot estimates of natural mortality (all size classes)
ggplot(dN.data, aes(x = Year2, y = Median, group = SizeClass, color = SizeClass)) + 
  geom_point() + geom_linerange(aes(x = Year2, ymin = lCI, ymax = uCI)) + 
  ylab('Natural mortality probability') + xlab('Year') + ylim(0,1) + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

## Plot estimates of natural mortality (large only)
ggplot(subset(dN.data, SizeClass == 'Large'), aes(x = Year, y = Median)) + 
  geom_point() + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI)) + 
  ylab('Natural mortality probability') + xlab('Year') + ylim(0,1) + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

## Plot estimates of hunting mortality (all size classes)
ggplot(dH.data, aes(x = Year2, y = Median, group = SizeClass, color = SizeClass)) + 
  geom_point() + geom_linerange(aes(x = Year2, ymin = lCI, ymax = uCI)) + 
  ylab('Hunting mortality probability') + xlab('Year') + 
  ylim(0,1) + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')


# 3) Postnatal survival over time #
#---------------------------------#

## Arrange data
s0.data <- as.data.frame(s0)
s0.data$Year <- years
rownames(s0.data) <- NULL

## Plot estimates of postnatal survival
ggplot(s0.data, aes(x = Year, y = Median)) + 
  geom_point(color = '#D9AD5B') + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI), color = '#D9AD5B') + 
  ylab('Postnatal survival probability') + xlab('Year') + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')


# 4) Breeding probability over time #
#-----------------------------------#

## Arrange data
pB.data <- as.data.frame(rbind(t(pB[,,1]), t(pB[,,2]), t(pB[,,3])))
pB.data$Year <- rep(years, 3)
pB.data$Year2 <- c(years-0.2, years, years+0.2)
pB.data$SizeClass <- rep(c('Small', 'Medium', 'Large'), each = length(years))
rownames(pB.data) <- NULL

## Plot estimates of breeding probability (all size classes)
ggplot(pB.data, aes(x = Year2, y = Median, group = SizeClass, color = SizeClass)) + 
  geom_point() + geom_linerange(aes(x = Year2, ymin = lCI, ymax = uCI)) + 
  ylab('Breeding probability') + xlab('Year') +  
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

## Plot estimates of breeding probability (large only)
ggplot(subset(pB.data, SizeClass == 'Large'), aes(x = Year, y = Median)) + 
  geom_point() + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI)) + 
  ylab('Bredding probability') + xlab('Year') +  
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')#, axis.text.x = element_text(angle = 45, vjust = 0.5))


# 5) Litter size over time #
#--------------------------#

## Arrange data
nF.data <- as.data.frame(rbind(t(nF[,,1]), t(nF[,,2]), t(nF[,,3])))
nF.data$Year <- rep(years, 3)
nF.data$Year2 <- c(years-0.2, years, years+0.2)
nF.data$SizeClass <- rep(c('Small', 'Medium', 'Large'), each = length(years))
rownames(nF.data) <- NULL

## Plot estimates of litter size (all size classes)
ggplot(nF.data, aes(x = Year2, y = Median, group = SizeClass, color = SizeClass)) + 
  geom_point() + geom_linerange(aes(x = Year2, ymin = lCI, ymax = uCI)) + 
  ylab('Litter size') + xlab('Year') +  
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')


# 6) Size class transitions over time #
#-------------------------------------#

## Arrange data
gSM.data <- as.data.frame(gSM)
gSM.data$Year <- years
rownames(gSM.data) <- NULL

gSL.data <- as.data.frame(gSL)
gSL.data$Year <- years
rownames(gSL.data) <- NULL

gML.data <- as.data.frame(gML)
gML.data$Year <- years
rownames(gML.data) <- NULL

## Plot estimates of growth parameters
plot_gSM <- ggplot(gSM.data, aes(x = Year, y = Median)) + 
  geom_point(color = '#586E2B') + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI), color = '#586E2B') + 
  ylab('Transition probability') + xlab('Year') + ylim(0,1) + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

plot_gSL <- ggplot(gSL.data, aes(x = Year, y = Median)) + 
  geom_point(color = '#586E2B') + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI), color = '#586E2B') + 
  ylab('Transition probability') + xlab('Year') + ylim(0,1) + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

plot_gML <- ggplot(gML.data, aes(x = Year, y = Median)) + 
  geom_point(color = '#586E2B') + geom_linerange(aes(x = Year, ymin = lCI, ymax = uCI), color = '#586E2B') + 
  ylab('Transition probability') + xlab('Year') + ylim(0,1) + 
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(linetype = 'dotted'), legend.position = 'none')

gridExtra::grid.arrange(plot_gSM, plot_gSL, plot_gML, ncol = 3)
 

# 7) Population structure over time #
#-----------------------------------#

## Arrange data
Prop.data <- as.data.frame(rbind(marProp, octProp, HProp))
Prop.data$Year <- rep(years, 3)
Prop.data$Period <- rep(c('March census', 'October census', 'Harvest'), each = length(years))
rownames(Prop.data) <- NULL

Prop.data <- reshape2::melt(Prop.data, id.vars = c('Year', 'Period'))
colnames(Prop.data) <- c('Year', 'Period', 'SizeClass', 'Proportion')

Prop.data$Period <- factor(Prop.data$Period, levels = c('March census', 'October census', 'Harvest'))

## Plot estimates of population structure
ggplot(Prop.data, aes(x = Year, y = Proportion, group = SizeClass)) +
  geom_bar(aes(fill = SizeClass), stat = 'identity') + 
  facet_wrap(~Period, ncol = 1) + 
  scale_fill_manual(values = c('grey60', 'grey82', 'black')) +
  scale_x_continuous(breaks = seq(1991, 2016, by = 5)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = '#D9AD5B'))


# 8) Natural and hunting mortality - densities #
#----------------------------------------------#

## Arrange data
dN.ddata <- data.frame(Estimate_m = c(out.mat[,'Mu.mN[1]'],
                                      out.mat[,'Mu.mN[2]'],
                                      out.mat[,'Mu.mN[3]']),
                       Estimate_d = c(1-exp(-out.mat[,'Mu.mN[1]']),
                                    1-exp(-out.mat[,'Mu.mN[2]']),
                                    1-exp(-out.mat[,'Mu.mN[3]'])),
                       SizeClass = rep(c('Small','Medium','Large'), each = nrow(out.mat)))

dH.ddata <- data.frame(Estimate_m = c(out.mat[,'Mu.mH[1]'],
                                      out.mat[,'Mu.mH[2]'],
                                      out.mat[,'Mu.mH[3]']),
                       Estimate_d = c(1-exp(-out.mat[,'Mu.mH[1]']),
                                      1-exp(-out.mat[,'Mu.mH[2]']),
                                      1-exp(-out.mat[,'Mu.mH[3]'])),
                       SizeClass = rep(c('Small','Medium','Large'), each = nrow(out.mat)))

## Plot posterior distributions

# Natural mortality hazard rate (all size classes)
ggplot(dN.ddata, aes(x = Estimate_m, group = SizeClass, color = SizeClass)) + 
  geom_density(fill = NA) + 
  ylab('Density') + xlab('Estimate') + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')

# Natural mortality probability (all size classes)
ggplot(dN.ddata, aes(x = Estimate_d, group = SizeClass, color = SizeClass)) + 
  geom_density(fill = NA, size = 0.75) + 
  ylab('Density') + xlab('Estimate') + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')

# Hunting mortality hazard rate (all size classes)
ggplot(dH.ddata, aes(x = Estimate_m, group = SizeClass, color = SizeClass)) + 
  geom_density(fill = NA) + 
  ylab('Density') + xlab('Estimate') + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')

# Hunting mortality probability (all size classes)
ggplot(dH.ddata, aes(x = Estimate_d, group = SizeClass, color = SizeClass)) + 
  geom_density(fill = NA) + 
  ylab('Density') + xlab('Estimate') + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')


# 9) Breeding probability - densities #
#-------------------------------------#

## Arrange data
pB.ddata <- data.frame(
  Estimate = c(out.mat[,'Mu.pB[1, 1]'],
               out.mat[,'Mu.pB[1, 2]'],
               out.mat[,'Mu.pB[1, 3]'],
               out.mat[,'Mu.pB[2, 1]'],
               out.mat[,'Mu.pB[2, 2]'],
               out.mat[,'Mu.pB[2, 3]'],
               out.mat[,'Mu.pB[3, 1]'],
               out.mat[,'Mu.pB[3, 2]'],
               out.mat[,'Mu.pB[3, 3]']),
  SizeClass = rep(c('Small', 'Medium', 'Large'), each = nrow(out.mat)*3),
  AcornAbundance = rep(rep(c('N', 'AA', 'A'), each = nrow(out.mat)), 3)
)


## Plot posterior distributions for breeding probability (all size classes and acorn abundances)
ggplot(pB.ddata) + geom_density(aes(x = Estimate, color = AcornAbundance), fill = NA, alpha = 0.5) + 
  scale_fill_manual(values = c('#44A694', 'black', '#D9AD5B')) + 
  scale_color_manual(values = c('#44A694', 'black', '#D9AD5B')) + 
  ylab('Density') + xlab('Estimate') + 
  facet_wrap(~SizeClass, scale = 'free_y') + theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none')


# 10) Recapture and recovery probability - densities #
#----------------------------------------------------#

## Arrange data
pp.ddata <- data.frame(Estimate = c(out.mat[,'Mu.pp[1]'],
                                    out.mat[,'Mu.pp[2]'],
                                    out.mat[,'Mu.pp[3]']),
                       SizeClass = rep(c('Small','Medium','Large'), each = nrow(out.mat)))

ll.ddata <- data.frame(Estimate = out.mat[,'Mu.ll'])

## Plot posterior distributions

# Recapture rate (all size classes)
ggplot(pp.ddata, aes(x = Estimate, group = SizeClass, color = SizeClass)) + 
  geom_density(fill = NA) + 
  ylab('Density') + xlab('Estimate') + 
  scale_color_manual(values = c('black', 'grey82', 'grey60')) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')

# Recovery rate (size-class-independent)
ggplot(ll.ddata, aes(x = Estimate)) + 
  geom_density(fill = NA) + 
  ylab('Density') + xlab('Estimate') + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = 'none')

