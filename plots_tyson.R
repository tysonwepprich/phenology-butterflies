#Plot stopover model results
library(msm)


setwd('C:/Users/Tyson/Dropbox/SESYNC ZGroup Summer2015/Ohio')

load(file = "fit_AZ_3M.RData")

count_array <- readRDS('count_array.rds')
cov_array <- readRDS('covariates_array.rds')
species_list <- readRDS('top20species.rds')

# source('FunctionsFixedForUnivoltineCase.R')
source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

counts <- count_array[,,4] #select species here, number corresponds to row in species_list

# select sites with enough individuals counted
siteRows <- which(rowSums(counts, na.rm = TRUE) >= 10)
counts <- counts[siteRows, ]
counts[is.na(counts)] <- -1


# Covariate (latitude)
lat <- rowMeans(cov_array[,,5], na.rm = TRUE)
lat <- lat[siteRows]
cov.w <- cov.mu <- scale(lat)[,1]


#---  for selected model ---#
FittedVal <- FittedVal.f(temp)
# Residual deviance
ResDev <- ResDev.f(FittedVal)
# DOF
DOF <- length(c(counts[counts > -1])) - temp$npar
# Residuals

overdis <- ResDev/DOF

VC <- solve(temp$Hessian)

#Eleni's plots
#plot N estimates with 95% CI

par(mar=c(5.1, 4.4, 4.1, 5.1))
library(plotrix)
plotCI((temp$N.est),ui=temp$N.up,li=temp$N.low,ylab="N (Np)",xlab="Site",cex.axis=1.1,cex.lab=1.4,axes=F)
axis(1,seq(1,S,1),tcl=0.5)
axis(2,seq(0,max(temp$N.up),50),tcl=0.5)
box()       

temp$phi.low[1,1,1]
temp$phi.est[1,1,1]
temp$phi.up[1,1,1]
# 
# plot(seq(1, T, 1), seq(0, 1, length = T),  "n", ylab = "Emegence curve", xlab = "Time")
# for(j in order(cov.mu)) lines(temp$betta.est[j,], col = j)


#plot counts vs fitted values for each site
# 
# for (i in order(cov.mu)){
#   plot(counts[i,], type = "l")
#   lines(FittedVal[i,], lty = 2)
# }

#Tyson's plots
S <- dim(counts)[1]
t <- dim(counts)[2]
df_all <- data.frame()
for (i in 1:S){
  df <- data.frame(Week = 1:t, Site = rep(i, t))
  df$Fit <- FittedVal[i,]
  df$Count <- counts[i,]
  df$Count[df$Count == -1] <- NA
  df$Phen <- temp$betta.est[i,]
  df$Latitude <- -cov.mu[i]
  
  df_all <- rbind(df_all, df)
}

library(dplyr)
df_all <- df_all %>% group_by(Site) %>% mutate(Tot = sum(Count, na.rm = TRUE)) %>% filter(Tot > 10)

library(ggplot2)

#plot fitted phenology by site, colored by latitude
c <- ggplot(data = df_all, aes(x = Week, y = Phen, group = Site)) + geom_line(aes(color = Latitude), size = .8) 
c + theme_bw() +
  scale_colour_gradient2(name = "Latitude", midpoint = mean(range(df_all$Latitude)), low = "red", mid = "yellow", high = "blue") 

#plot fitted values and actual counts for all sites
#Scaled latitude is multiplied by -1, so Southern sites at bottom of the plot
d <- ggplot(data = df_all, aes(x = Week, y = Count, group = Site)) + geom_point() + 
  facet_wrap(Latitude ~ Site, ncol = 4, scales = "free_y") + geom_line(aes(x = Week, y = Fit)) + theme_bw()
d




















