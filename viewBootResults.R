#Visualize results from bootstrapped stopover models
# setwd('R/StopoverModel')
setwd('C:/Users/Tyson/Dropbox/SESYNC ZGroup Summer2015/Ohio')


source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

#load data
bootFits <- readRDS("LWS_Boot_incomp.rds")
bootIdxs <- readRDS("LWS_id_incomp.rds")
timeList <- readRDS("LWStime.rds")

#fitting times for each bootstrap
times <- unlist(lapply(timeList, as.double, units = "mins"))



nBoots <- length(bootFits)
nSites <- length(bootFits[[1]]$N.est)
nWeeks <- dim(bootFits[[1]]$betta.est)[2]
M <- dim(bootFits[[1]]$mu.est)[2]
qp <- length(bootFits[[1]]$cest)

# I copied the following from Rebecca's code so you might want to change it, especially for extracting the estimates for the coefficients for p (which can be more than 2 in your case)
bootB0 = array(NA,c(nBoots, M)) #mean arrival times
bootB1 = array(NA,c(nBoots,1))
bootC.est = array(NA,c(nBoots, qp)) # here for example you want to extract qp+1 sets of estimates
# bootC1 = array(NA,c(nBoots,1))
bootD0 = array(NA, c(nBoots, length(plg(rep(1/M, M))))) #unsure if this length is correct, line 115 of function code
bootD1 = array(NA, c(nBoots, 1))
bootSigma = array(NA,c(nBoots, M)) #would only need 1 if sigma is 'hom'
bootPhi = array(NA,c(nBoots,1))
bootNs = array(NA,c(nBoots, nSites))
bootMus = array(NA,c(nBoots, nSites, M))
bootWs = array(NA, c(nBoots, nSites, M))
bootBetta = array(NA, c(nBoots, nSites, nWeeks))
bootLL = array(NA,c(nBoots,1))
for (b in 1:nBoots) {
  bootLL[b] = bootFits[[b]]$ll.val
  bootB0[b, ] = bootFits[[b]]$b0.est
  bootB1[b] = bootFits[[b]]$b1.est
  bootC.est[b, ] = bootFits[[b]]$cest
  bootD0[b, ] = bootFits[[b]]$d0.est
  bootD1[b] = bootFits[[b]]$d1.est
  bootSigma[b, ] = bootFits[[b]]$sigma.est[1, ]
  bootPhi[b] = bootFits[[b]]$phi.est[1,1,1]
  for (s in 1:length(bootFits[[1]]$N.est)) {
    idxsThisSite = which(bootIdxs[b,]==s)
    if (length(idxsThisSite)>0) {
      bootNs[b,s] = bootFits[[b]]$N.est[idxsThisSite[1]]
      bootMus[b,s, ] = bootFits[[b]]$mu.est[idxsThisSite[1], ]
      bootWs[b,s, ] = bootFits[[b]]$w.est[idxsThisSite[1], ]
      bootBetta[b,s, ] = bootFits[[b]]$betta.est[idxsThisSite[1], ]
    
    }
  }
}


#plot N estimates for each site
boxplot(x = as.list(as.data.frame(bootNs)), pars = list(ylim = c(0, 1000)))

#plot mu estimates for each site

#interesting, seems like it doesn't order the mu's by time, so 1st, 2nd, 3rd 
#generations get interchanged
boxplot(x = as.list(as.data.frame(bootMus[,,1])))
boxplot(x = as.list(as.data.frame(bootMus[,,2])))
#happens in about half the bootstraps
#mu for second brood is earlier than that of first brood
which(rowMeans(bootMus[,,1], na.rm = TRUE) > rowMeans(bootMus[,,2], na.rm = TRUE))

#however, histogram of all mus seems ok
hist(bootMus[,,1][bootMus[,,1]<50], breaks = 100)
hist(bootMus[,,2][bootMus[,,2]<50], breaks = 100)

#by site
hist(bootMus[,1,], breaks = 50)

#weights for each generation
hist(rowMeans(bootWs[,,1], na.rm = TRUE), breaks = 25)
hist(rowMeans(bootWs[,,2], na.rm = TRUE), breaks = 25)

#rearrange Mus and Ws so that broods correctly separated post-fitting
bootMus_adj = array(NA,c(nBoots, nSites, M))
bootWs_adj = array(NA, c(nBoots, nSites, M))
bootSigma_adj = array(NA,c(nBoots, M)) #would only need 1 if sigma is 'hom'
for (i in 1:nBoots){
  colMus <- colMeans(bootMus[i,,], na.rm = TRUE)
  sortMus <- sort(colMus)
  for (m in 1:M){
    bootMus_adj[i,,m] <- bootMus[i,, which(colMus == sortMus[m])]
    bootWs_adj[i,,m] <- bootWs[i,, which(colMus == sortMus[m])]
    bootSigma_adj[i,m] <- bootSigma[i,which(colMus == sortMus[m])]
  }
}

#when adjusted, mostly OK, except for a few where 1st generation way too early
hist(rowMeans(bootMus_adj[,,1], na.rm = TRUE), breaks = 50)
hist(rowMeans(bootMus_adj[,,2], na.rm = TRUE), breaks = 50)

#weights averaged by site before and after rearranging Mus and Ws
#adjusted shows that first brood more populous than second.
colMeans(bootWs[,,1], na.rm = TRUE)
colMeans(bootWs[,,2], na.rm = TRUE)
#versus
colMeans(bootWs_adj[,,1], na.rm = TRUE)
colMeans(bootWs_adj[,,2], na.rm = TRUE)

hist(rowMeans(bootWs_adj[,,1], na.rm = TRUE), breaks = 25)
hist(rowMeans(bootWs_adj[,,2], na.rm = TRUE), breaks = 25)


#however, I think switching the broods makes the covariate estimates off
#this removes two big outliers
#then plots Mu estimates for each brood against D0 estimates
#looks like switching brood order changes the sign of D0
plot(rowMeans(bootMus[-c(14,35),,1], na.rm = TRUE), bootD0[-c(14,35)])
plot(rowMeans(bootMus[-c(14,35),,2], na.rm = TRUE), bootD0[-c(14,35)])




#histogram of other estimates

hist(bootSigma_adj[-14,1], breaks = 25)
hist(bootSigma_adj[-14,2], breaks = 25)


#loglik
hist(bootLL, breaks = 25)

#mean arrival times with covariate
hist(bootB0[,1], breaks = 25)
hist(bootB0[,2], breaks = 25)

hist(bootB1, breaks = 25)

#detection probability with covariates
hist(bootC.est[,1], breaks = 25)
hist(bootC.est[,2], breaks = 25)
hist(bootC.est[,3], breaks = 25)
colMeans(bootC.est)

#brood weights with covariate
hist(bootD0[abs(bootD0)<5], breaks = 25)
hist(bootD1[abs(bootD1)<5], breaks = 25)

hist(bootSigma[,1], breaks = 25)
hist(bootSigma[,2], breaks = 25)

hist(bootPhi, breaks = 25)


