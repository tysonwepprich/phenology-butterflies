#Code stopped after first successful convergence: 
#Error in mLLMixtCounts.fit(p.m, w.m, mu.m, sigma.m, phi.m) : 
#object 'invect' not found
#replaced it with 'outvect', seems to fix problem



rm(list=ls(all=TRUE))

setwd('C:/Users/Tyson/Dropbox/SESYNC ZGroup Summer2015/Ohio')

count_array <- readRDS('count_array.rds')
cov_array <- readRDS('covariates_array.rds')
species_list <- readRDS('top20species.rds')




# source('FunctionsFixedForUnivoltineCase.R')
source('FunctionsFixedForUnivoltineCaseMultipleDetectionCovariates.R')

counts <- count_array[,,12] #select species here, number corresponds to row in species_list

# select sites with enough individuals counted
siteRows <- which(rowSums(counts, na.rm = TRUE) >= 10)
counts <- counts[siteRows, ]
counts[is.na(counts)] <- -1

M <- 2
S <- dim(counts)[1]
K <- T <- dim(counts)[2]

# Covariate weather for p
# cov.p <- matrix(data_samp$TEMP,nrow=S,ncol=K,byrow=TRUE)
#cov.p now needs to be a S by T by qp (number of covariates) array
covs <- c(2,4) #Select detection covariate here (1:5 possible)
if (length(covs) > 1) cov.p <- cov_array[,,covs]
if (length(covs) == 1) cov.p <- array(data = cov_array[,,covs], dim = c(dim(cov_array[,,covs]), 1))
cov.p <- cov.p[siteRows, , , drop = FALSE]
# qp is the number of covariates for p
qp <- dim(cov.p)[3]
for(q in 1:qp) cov.p[,,q] <- scale(cov.p[,,q])[1:S,1:K]
cov.p[is.na(cov.p)] <- -1

#Aside, ignore for now
#What about a covariate of temperature ^ 2 for detection?
# library(abind)
# abind(la[[1]],la[[2]],la[[3]],along = 1)



# Time covariate for phi
cov.phi <-  matrix(((1:(K-1))-mean(1:(K-1)))/sqrt(var(1:(K-1))),S,K-1,byrow=TRUE) 

# Covariate (latitude)
lat <- rowMeans(cov_array[,,5], na.rm = TRUE)
lat <- lat[siteRows]
cov.w <- cov.mu <- scale(lat)[,1]


# add bootstrap

transCounts = counts
trans.cov.phi = cov.phi
trans.cov.w = cov.w
trans.cov.mu = cov.mu
trans.cov.p = cov.p


nBoots = 100
# bootFits = list()
# bootIdxs = array(0,c(nBoots,dim(counts)[1]))
# bootTime <- list()

bootFits <- readRDS("GSFBoot.rds")
bootIdxs <- readRDS("GSFIDxs.rds")
bootTime <- readRDS("GSFTime.rds")


# 8 bootstrap replicates took 2 days!
for (b in 48:nBoots) {
  
  print(b)
  startTime <- Sys.time()
  
  bSites = sample(dim(counts)[1],replace=T)
  bootIdxs[b,] = bSites
  counts = transCounts[bSites,]
  cov.w = trans.cov.w[bSites]
  cov.mu = trans.cov.mu[bSites]
  cov.p = trans.cov.p[bSites,,, drop = FALSE]

## MODEL RUNNING ######################################

########################################################
p.m <- "cov"
w.m <- "cov"
mu.m <- "cov"
sigma.m <- "het"
# phi.m <- "quadr.t"
phi.m <- "const"

#Think about putting it some try loop, so fits ending in errors don't stop program
#Maybe stop while loop when 3-5 temp.ll are not NA to choose max
# TryFit <- function(){
#   res <- lapply(1:10, function(i) try(mLLMixtCounts.fit(p.m,w.m,mu.m,sigma.m,phi.m), TRUE))
# }



temp.fit <- list()
temp.ll <- rep(NA,5)
k <- 1
check <- 1
while(k<6){
  start.list <- start.fun(p.m,w.m,mu.m,sigma.m,phi.m)
  #start.list$d0 <- 1
  pars.start <- c(start.list$N,start.list$cvec, start.list$d0, start.list$d1, start.list$b0, start.list$b1,start.list$sigma,  start.list$a0, start.list$a1, start.list$a2) #this line remains the same for all models
  st1 <- Sys.time()
  try(temp.fit[[k]] <- mLLMixtCounts.fit(p.m,w.m,mu.m,sigma.m,phi.m), TRUE)
  et1 <- Sys.time()
  if (is.na(temp.fit[[k]]$ll.val) == FALSE){
    temp.ll[k] <- temp.fit[[k]]$ll.val
    k <- k+1
  }else{
    check <- check + 1
  }
  if (check > 5) next
}

tempchoose <- min(c(1:5)[which(temp.ll==max(temp.ll,na.rm=T))])

# temp <-temp.fit[[tempchoose]] 

bootFits[[b]] <- temp.fit[[tempchoose]]
bootTime[[b]] <- Sys.time() - startTime

saveRDS(bootFits, "GSFBoot.rds")
saveRDS(bootIdxs, "GSFIDxs.rds")
saveRDS(bootTime, "GSFTime.rds")
} # b



# I copied the following from Rebecca's code so you might want to change it, especially for extracting the estimates for the coefficients for p (which can be more than 2 in your case)
bootB0 = array(NA,c(nBoots, M)) #mean arrival times
bootB1 = array(NA,c(nBoots,1))
bootC.est = array(NA,c(nBoots, qp + 1)) # here for example you want to extract qp+1 sets of estimates
# bootC1 = array(NA,c(nBoots,1))
bootD0 = array(NA, c(nBoots, length(plg(rep(1/M, M))))) #unsure if this length is correct, line 115 of function code
bootD1 = array(NA, c(nBoots, 1))
bootSigma = array(NA,c(nBoots, M)) #would only need 1 if sigma is 'hom'
bootPhi = array(NA,c(nBoots,1))
bootNs = array(NA,c(nBoots, dim(counts)[1]))
bootMus = array(NA,c(nBoots, dim(counts)[1], M))
bootWs = array(NA, c(nBoots, dim(counts)[1], M))
bootBetta = array(NA, c(nBoots, dim(counts)[1], dim(counts)[2]))
for (b in 1:nBoots) {
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
      
#       if (s==1) {
#         bootBetaSite1[b,] = bootFits[[b]]$betta.est[idxsThisSite[1],]
#       }
    }
  }
}


# save(temp,file=paste("C://Users/ed234/Dropbox/so_mix/weeks/latestCI/fit_",i_species,"_2010_s",S,"_12_SO_temp.RData",sep=""))
save(temp, file = "fit_AZ_3M.RData")

###############################################################


# load(file=paste("C://Users/ed234/Dropbox/so_mix/weeks/latestCI/fit_",i_species,"_2010_s",S,"_6_SO_temp.RData",sep=""))
load(file = "fit_LWS.RData")
cat(round(temp$ll.val,2),"&",temp$npar,"&",round(2*temp$npar-2*temp$ll.val,2),"\\\\\n")	

# Best model
fit3 <- temp 

# Stopover duration
duration <- Stopover.f(temp)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,duration[order(cov.w)],type="l",xlab="Northing (km)",ylab="Stopover duration",cex.axis=1.1,cex.lab=1.4)


######################################################### RESULTS #####################################################
#---  for selected model ---#
FittedVal <- FittedVal.f(temp)
# Residual deviance
ResDev <- ResDev.f(FittedVal)
# DOF
DOF <- length(c(counts[counts > -1]))-fit3$npar
# Residuals

overdis <- ResDev/DOF

library(msm)

VC <- solve(fit3$Hessian)


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

plot(seq(1, T, 1), seq(0, 1, length = T),  "n", ylab = "Emegence curve", xlab = "Time")
for(j in order(cov.mu)) lines(temp$betta.est[j,], col = j)



temp_temp <- seq(50, 95, 0.5)
temp_temp_sc <- scale(temp_temp)[,1]
plot(temp_temp, expo(temp$c0.est + temp$c1.est*temp_temp_sc), type = "l", xlab = "Temperature", ylab = "Detection probability")
plot(cov.p[is.na(cov.p)==F], temp$p.est[is.na(cov.p)==F])

plot(sort(cov.p[1,][is.na(cov.p[1,])==F]), c(temp$p.est[1,][is.na(cov.p[1,])==F])[order(cov.p[1,][is.na(cov.p[1,])==F])])

for(i in 1:S){
  for(t in 1:T){
    if(is.na(cov.p[i,t])==F) points(cov.p[i,t], temp$p.est[i,t])
  }
}
  

#--- Estimates and standard errors ---#

Nlow <- exp(log(fit3$N.est)-1.96*sqrt(diag(VC)[1:S]*overdis));Nup <- exp(log(fit3$N.est)+1.96*sqrt(diag(VC)[1:S]*overdis))


postscript("P:/Documents/MixtCounts/Paper/Figures/Ncomp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
library(plotrix)
plotCI((fit3$N.est),ui=Nup,li=Nlow,ylab="N (Np)",xlab="Site",cex.axis=1.1,cex.lab=1.4,axes=F)
axis(1,seq(1,S,1),tcl=0.5)
axis(2,seq(0,max(Nup),50),tcl=0.5)
points(fit1$N.est, pch=19)
dev.off()

pdf("D:/matechou/Documents/Talks/NCSE2013/Ncomp.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
library(plotrix)
plotCI((fit3$N.est),ui=Nup,li=Nlow,ylab="N",xlab="Site",cex.axis=1.1,cex.lab=1.4,axes=F)
axis(1,seq(1,S,1),tcl=0.5)
axis(2,seq(0,max(Nup),50),tcl=0.5)
points(fit1$N.est, pch=19)
dev.off()
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for p ---#

round(c(fit3$c0.est,sqrt(diag(VC)[S+1]*overdis)),3)

round(c(fit3$c1.est,sqrt(diag(VC)[S+2]*overdis)),3)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for w ---#

round(c(fit3$d0.est,sqrt(diag(VC)[S+3]*overdis)),3)

round(c(fit3$d1.est,sqrt(diag(VC)[S+4]*overdis)),3)

w1low<-c();w1up <- c()

for(i in 1:S){
		form <- sprintf("~x1 + x2*%f",cov.w[i])
		w1low[i] <- expo(fit3$d0.est + fit3$d1.est * cov.w[i] - 1.96*deltamethod(as.formula(form),c(fit3$d0.est,fit3$d1.est),VC[(53):(54),(53):(54)]*overdis))
		w1up[i] <- expo(fit3$d0.est + fit3$d1.est * cov.w[i] + 1.96*deltamethod(as.formula(form),c(fit3$d0.est,fit3$d1.est),VC[(53):(54),(53):(54)]*overdis))
				}


postscript("P:/Documents/MixtCounts/Paper/Figures/wlat_temp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit3$w.est[,1][order(cov.w)],type="l",ylim=c(0,1),xlab="Northing (km)",
ylab=expression(w[1]),cex.axis=1.1,cex.lab=1.4,axes=F,xlim=c(0,1000))
axis(1,seq(0,1000,100),tcl=0.5)
axis(2,seq(0,1,0.1),tcl=0.5)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1low[order(cov.w)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1up[order(cov.w)],lty=2,type="l")
dev.off()


pdf("D:/matechou/Documents/Talks/NCSE2013/wlat_temp.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit3$w.est[,1][order(cov.w)],type="l",ylim=c(0,1),xlab="Northing (km)",
ylab=expression(w[1]),cex.axis=1.1,cex.lab=1.4,axes=F,xlim=c(0,1000))
axis(1,seq(0,1000,100),tcl=0.5)
axis(2,seq(0,1,0.1),tcl=0.5)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1low[order(cov.w)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1up[order(cov.w)],lty=2,type="l")
dev.off()


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for mu ---#
round(c(fit3$b0.est[1],sqrt(diag(VC)[55]*overdis)),3)

round(c(fit3$b0.est[2],sqrt(diag(VC)[56]*overdis)),3)

round(c(fit3$b1.est,sqrt(diag(VC)[57]*overdis)),3)


mulow<-matrix(NA,S,M); muup<-matrix(NA,S,M)

		for(i in 1:S)
			{
			for(m in 1:M)
				{
				form <- sprintf("~x1 + x2*%f",cov.mu[i])
				mulow[i,m] <- exp(fit3$b0.est[m] + fit3$b1.est * cov.mu[i] - 1.96*deltamethod(as.formula(form),c(fit3$b0.est[m],fit3$b1.est),VC[c((54+m),(57)),c((54+m),(57))]*overdis))
				muup[i,m] <- exp(fit3$b0.est[m] + fit3$b1.est * cov.mu[i] + 1.96*deltamethod(as.formula(form),c(fit3$b0.est[m],fit3$b1.est),VC[c((54+m),(57)),c((54+m),(57))]*overdis))
				}
             	}


postscript("P:/Documents/MixtCounts/Paper/Figures/mulat_temp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))

plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit3$mu.est[,1][order(cov.mu)],type="l",ylim=c(1,29),xlab="Northing (km)",
ylab=expression(mu),cex.axis=1.1,cex.lab=1.4,axes=F,xlim=c(0,1000))
axis(1,seq(0,1000,100),tcl=0.5)
axis(2,seq(0,30,1),tcl=0.5)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow[,1][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup[,1][order(cov.mu)],lty=2,type="l")

lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit3$mu.est[,2][order(cov.mu)],type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow[,2][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup[,2][order(cov.mu)],lty=2,type="l")

dev.off()

pdf("D:/matechou/Documents/Talks/NCSE2013/mulat_temp.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))

plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit3$mu.est[,1][order(cov.mu)],type="l",ylim=c(1,29),xlab="Northing (km)",
ylab=expression(mu),cex.axis=1.1,cex.lab=1.4,axes=F,xlim=c(0,1000))
axis(1,seq(0,1000,100),tcl=0.5)
axis(2,seq(0,30,1),tcl=0.5)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow[,1][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup[,1][order(cov.mu)],lty=2,type="l")

lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit3$mu.est[,2][order(cov.mu)],type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow[,2][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup[,2][order(cov.mu)],lty=2,type="l")

dev.off()


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- sigma ---#

round(c(fit3$sigma.est[1,1],sqrt(diag(VC)[58]*overdis)),3)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for phi ---#

round(c(fit3$inter.est,sqrt(diag(VC)[59]*overdis)),3)

round(c(fit3$slope.est,sqrt(diag(VC)[60]*overdis)),3)

round(c(fit3$slope2.est,sqrt(diag(VC)[61]*overdis)),3)

philow <- array(0,c(T-1,T-1,S)); phiup <- array(0,c(T-1,T-1,S))

for(i in 1:S){	
		for(j in 1:(T-1)){
			form <- sprintf("~x1 + x2*%f + x3*%g",cov.phi[i,j],cov.phi[i,j]^2)
			philow[,j,i] <- expo(c(fit3$inter.est + fit3$slope.est*cov.phi[i,j] + fit3$slope2.est*cov.phi[i,j]^2) - 1.96*deltamethod(as.formula(form),c(fit3$inter.est,fit3$slope.est,fit3$slope2.est),VC[59:61,59:61]*overdis))
			phiup[,j,i] <- expo(c(fit3$inter.est + fit3$slope.est*cov.phi[i,j] + fit3$slope2.est*cov.phi[i,j]^2) + 1.96*deltamethod(as.formula(form),c(fit3$inter.est,fit3$slope.est,fit3$slope2.est),VC[59:61,59:61]*overdis))
					  }
			philow[,,i][lower.tri(philow[,,i])]<-0; phiup[,,i][lower.tri(phiup[,,i])]<-0;
			}

postscript("P:/Documents/MixtCounts/Paper/Figures/phitime.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(fit3$phi.est[1,,1],type="l",xlab="Week",ylab=expression(phi),ylim=c(0,1),cex.axis=1.1,cex.lab=1.4,axes=F)
axis(1,seq(1,K-1,1),tcl=0.5)
axis(2,seq(0,1,0.1),tcl=0.5)
lines(philow[1,,1],type="l",lty=2)
lines(phiup[1,,1],type="l",lty=2)
dev.off()


pdf("D:/matechou/Documents/Talks/NCSE2013/phitime.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(fit3$phi.est[1,,1],type="l",xlab="Week",ylab=expression(phi),ylim=c(0,1),cex.axis=1.1,cex.lab=1.4,axes=F)
axis(1,seq(1,K-1,1),tcl=0.5)
axis(2,seq(0,1,0.1),tcl=0.5)
lines(philow[1,,1],type="l",lty=2)
lines(phiup[1,,1],type="l",lty=2)
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- observed together with fitted counts ---#

countsNA <- counts
countsNA[countsNA == -1] <- NA
postscript("P:/Documents/MixtCounts/Paper/Figures/countsOE1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 11, width = 8)
par(mfrow=c(5,5),oma=c(.5,.5,.5,.5),mar=c(4,4,.5,.5))		
for(i in 1:25){
	plot(countsNA[i,],type="o",xlim=c(1,K),xlab="",ylab="",cex.axis=1.1)
	lines(1:K,FittedVal[i,],lty=2)
	mtext("Week", side=1, line=2.5,cex=.8)
	mtext("Count", side=2, line=2.5,cex=.8)
	}
dev.off()
postscript("P:/Documents/MixtCounts/Paper/Figures/countsOE2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 11, width = 8)
par(mfrow=c(5,5),oma=c(.5,.5,.5,.5),mar=c(4,4,.5,.5))	
for(i in 26:50){
	plot(countsNA[i,],type="o",xlim=c(1,K),xlab="",ylab="",cex.axis=1.1,axes=F)
axis(1,seq(1,K,2))
axis(2,tcl=0.5)
	lines(1:K,FittedVal[i,],lty=2)
	mtext("Week", side=1, line=2.5,cex=.8)
	mtext("Count", side=2, line=2.5,cex=.8)
	}
dev.off()
#-----------


pdf("D:/matechou/Documents/Talks/NCSE2013/countsOE1.pdf", height = 8, width = 8)
par(mfrow=c(5,5),oma=c(.5,.5,.5,.5),mar=c(4,4,.5,.5))		
for(i in 1:25){
	plot(countsNA[i,],type="o",xlim=c(1,K),xlab="",ylab="",cex.axis=1.1)
	lines(1:K,FittedVal[i,],lty=2)
	mtext("Week", side=1, line=2.5,cex=.8)
	mtext("Count", side=2, line=2.5,cex=.8)
	}
dev.off()
pdf("D:/matechou/Documents/Talks/NCSE2013/countsOE2.pdf", height = 8, width = 8)
par(mfrow=c(5,5),oma=c(.5,.5,.5,.5),mar=c(4,4,.5,.5))	
for(i in 26:50){
	plot(countsNA[i,],type="o",xlim=c(1,K),xlab="",ylab="",cex.axis=1.1,axes=F)
axis(1,seq(1,K,2))
axis(2,tcl=0.5)
	lines(1:K,FittedVal[i,],lty=2)
	mtext("Week", side=1, line=2.5,cex=.8)
	mtext("Count", side=2, line=2.5,cex=.8)
	}
dev.off()

##----------------------------------------------------------------------------------------------------------------------------------------------------------
#--- estimated beta parameters ---#

# Plots with covariate

latitudes <- (seq(floor(min(data_samp$NORTH)/10000)*10000,ceiling(max(data_samp$NORTH)/10000)*10000,1000)- attr(cov.mu,"scaled:center"))/ attr(cov.mu,"scaled:scale")

betta.lat <-  matrix(0,nrow=length(latitudes),ncol=K)
mu.lat1 <- w.lat <- matrix(NA,nrow=length(latitudes), ncol=2)

for(i in 1:length(latitudes)){
	mu.lat1[i,] <- c(exp(fit3$b0.est[1] + fit3$b1.est*latitudes[i]),exp(fit3$b0.est[2] + fit3$b1.est*latitudes[i]))
	w.lat[i,] <- lgp(fit3$d0.est + fit3$d1.est*latitudes[i])
	for(m in 1:M){
		betta.lat[i,] <- betta.lat[i,]+c(w.lat[i,m]*c(pnorm(1,mean=mu.lat1[i,m],sd=fit3$sigma.est[1,1]),
		pnorm(2:(K-1),mean=mu.lat1[i,m],sd=fit3$sigma.est[1,1])-pnorm(1:(K-2),mean=mu.lat1[i,m],sd=fit3$sigma.est[1,1]),		
		1-pnorm(K-1,mean=mu.lat1[i,m],sd=fit3$sigma.est[1,1])))
		}
	}

lats <- seq(floor(min(data_samp$NORTH)/10000)*10000,ceiling(max(data_samp$NORTH)/10000)*10000,1000)
	
postscript("P:/Documents/MixtCounts/Paper/Figures/bettalat.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 11, width = 8)
par(mfrow=c(4,3),oma=c(1,1,1,1),mar=c(4,5,1,1))
for(i in which(lats %in% seq(75000,950000,75000)))
{
      if(i!=701&i!=776&i!=851){ 
	plot(betta.lat[i,],ylab=expression(beta),ylim=c(0,.35),type="o",cex.lab=1.1,cex.main=0.9,cex.axis=1.1,xlab=NA)
	legend("topright",paste(lats[i]/1000,"(km)"),bty="n",cex=1.4)
					}
      if(i==701|i==776|i==851){ 
	plot(betta.lat[i,],ylab=expression(beta),ylim=c(0,.35),type="o",cex.lab=1.1,cex.main=0.9,cex.axis=1.1,xlab="Week")
	legend("topright",paste(lats[i]/1000,"(km)"),bty="n",cex=1.4)
					}
	}
dev.off()


pdf("D:/matechou/Documents/Talks/NCSE2013/bettalat.pdf", height = 8, width = 8)
par(mfrow=c(4,3),oma=c(1,1,1,1),mar=c(5,5,3,1))
for(i in which(lats %in% seq(75000,950000,75000))){
	plot(betta.lat[i,],ylab=expression(beta),xlab="Week",main=paste("Northing",lats[i]/1000,"(km)"),ylim=c(0,.35),type="o",cex.lab=1.4,cex.axis=1.1,cex.main=1.4,cex.axis=1.1)
	}
dev.off()

#--------------------------------- Mean stopover duration ----------------------------------------------#

plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,Stopover.f(fit3)[order(cov.mu)],type="l",ylim=c(1,5),xlab="Northing (km)",
ylab="Estimated mean stopover duration",cex.axis=1.1,cex.lab=1.4,axes=F,xlim=c(0,1000))
axis(1,seq(0,1000,100),tcl=0.5)
axis(2,seq(0,5,1),tcl=0.5)













#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################

#---  for 2nd best model ---#
FittedVal2 <- FittedVal.f(fit1)
# Residual deviance
ResDev2 <- ResDev.f(FittedVal2)
# DOF
DOF2 <- length(c(counts[counts > -1]))-fit1$npar
# Residuals

overdis2 <- ResDev2/DOF2

library(msm)

VC2 <- solve(fit1$Hessian)

#--- Estimates and standard errors ---#

Nlow2 <- exp(log(fit1$N.est)-1.96*sqrt(diag(VC2)[1:S]*overdis2));Nup2 <- exp(log(fit1$N.est)+1.96*sqrt(diag(VC2)[1:S]*overdis2))


postscript("P:/Documents/MixtCounts/Paper/Figures/Ncomp2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
library(plotrix)
plotCI((fit1$N.est),ui=Nup2,li=Nlow2,ylab="N",xlab="Site",cex.axis=1.1,cex.lab=1.4,ylim=c(0,400))
points(fit3$N.est, pch=19)
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for w ---#

round(c(fit1$d0.est,sqrt(diag(VC2)[51]*overdis2)),3)

round(c(fit1$d1.est,sqrt(diag(VC2)[52]*overdis2)),3)

w1low2<-c();w1up2 <- c()

for(i in 1:S){
		form <- sprintf("~x1 + x2*%f",cov.w[i])
		w1low2[i] <- expo(fit1$d0.est + fit1$d1.est * cov.w[i] - 1.96*deltamethod(as.formula(form),c(fit1$d0.est,fit1$d1.est),VC2[(51):(52),(51):(52)]*overdis2))
		w1up2[i] <- expo(fit1$d0.est + fit1$d1.est * cov.w[i] + 1.96*deltamethod(as.formula(form),c(fit1$d0.est,fit1$d1.est),VC2[(51):(52),(51):(52)]*overdis2))
				}


postscript("P:/Documents/MixtCounts/Paper/Figures/wlat_temp2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit1$w.est[,1][order(cov.w)],type="l",ylim=c(0,1),xlab="Northing (km)",ylab=expression(w),cex.axis=1.1,cex.lab=1.4)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1low2[order(cov.w)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,w1up2[order(cov.w)],lty=2,type="l")
dev.off()


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for mu ---#
round(c(fit1$b0.est[1],sqrt(diag(VC2)[53]*overdis2)),3)

round(c(fit1$b0.est[2],sqrt(diag(VC2)[54]*overdis2)),3)

round(c(fit1$b1.est,sqrt(diag(VC2)[55]*overdis2)),3)


mulow2<-matrix(NA,S,M); muup2<-matrix(NA,S,M)

		for(i in 1:S)
			{
			for(m in 1:M)
				{
				form <- sprintf("~x1 + x2*%f",cov.mu[i])
				mulow2[i,m] <- exp(fit1$b0.est[m] + fit1$b1.est * cov.mu[i] - 1.96*deltamethod(as.formula(form),c(fit1$b0.est[m],fit1$b1.est),VC2[c((52+m),(55)),c((52+m),(55))]*overdis2))
				muup2[i,m] <- exp(fit1$b0.est[m] + fit1$b1.est * cov.mu[i] + 1.96*deltamethod(as.formula(form),c(fit1$b0.est[m],fit1$b1.est),VC2[c((52+m),(55)),c((52+m),(55))]*overdis2))
				}
             	}


postscript("P:/Documents/MixtCounts/Paper/Figures/mulat_temp2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))

plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit1$mu.est[,1][order(cov.mu)],type="l",ylim=c(1,29),xlab="Northing (km)",ylab=expression(mu),cex.axis=1.1,cex.lab=1.4)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow2[,1][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup2[,1][order(cov.mu)],lty=2,type="l")

lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,fit1$mu.est[,2][order(cov.mu)],type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,mulow2[,2][order(cov.mu)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.mu)]/1000,muup2[,2][order(cov.mu)],lty=2,type="l")

dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- sigma ---#

round(c(fit1$sigma.est[1,1],sqrt(diag(VC2)[56]*overdis2)),3)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--- coefficients for phi ---#

round(c(fit1$inter.est,sqrt(diag(VC2)[57]*overdis2)),3)

round(c(fit1$slope.est,sqrt(diag(VC2)[58]*overdis2)),3)

round(c(fit1$slope2.est,sqrt(diag(VC2)[59]*overdis2)),3)

philow2 <- array(0,c(T-1,T-1,S)); phiup2 <- array(0,c(T-1,T-1,S))

for(i in 1:S){	
		for(j in 1:(T-1)){
			form <- sprintf("~x1 + x2*%f + x3*%g",cov.phi[i,j],cov.phi[i,j]^2)
			philow2[,j,i] <- expo(c(fit1$inter.est + fit1$slope.est*cov.phi[i,j] + fit1$slope2.est*cov.phi[i,j]^2) - 1.96*deltamethod(as.formula(form),c(fit1$inter.est,fit1$slope.est,fit1$slope2.est),VC2[57:59,57:59]*overdis2))
			phiup2[,j,i] <- expo(c(fit1$inter.est + fit1$slope.est*cov.phi[i,j] + fit1$slope2.est*cov.phi[i,j]^2) + 1.96*deltamethod(as.formula(form),c(fit1$inter.est,fit1$slope.est,fit1$slope2.est),VC2[57:59,57:59]*overdis2))
					  }
			philow2[,,i][lower.tri(philow2[,,i])]<-0; phiup2[,,i][lower.tri(phiup2[,,i])]<-0;
			}

postscript("P:/Documents/MixtCounts/Paper/Figures/phitime2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(fit1$phi.est[1,,1],type="l",xlab="Week",ylab=expression(phi),ylim=c(0,1),cex.axis=1.1,cex.lab=1.4)
lines(philow2[1,,1],type="l",lty=2)
lines(phiup2[1,,1],type="l",lty=2)
dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##############################################################################################################################################
#Not used

































































c0.out <- NULL; c1.out <- NULL
p.low <- matrix(NA, S, T); p.up <- matrix(NA, S, T)

if(p.m == "common"){ N.out <- exp(outvect[1:S]); par.index <- S
				  p.out <- matrix(1, nrow = S, ncol = T)
				} else{N.out <- exp(outvect[1:S])
					  c0.out <- outvect[S+1]; c1.out <- outvect[S+2]; par.index <- S + 2
					  p.out <- matrix(expo(c0.out + c1.out*cov.p), byrow = F, nrow = S, ncol = T)
for(i in 1:S){
	for(j in 1:T){
		form <- sprintf("~x1 + x2*%f",cov.p[i,j])
          p.low[i,j] <- expo(c0.out + c1.out*cov.p[i,j] - 1.96*deltamethod(as.formula(form),c(c0.out,c1.out),VC[(S+1):(S+2),(S+1):(S+2)]))
          p.up[i,j] <- expo(c0.out + c1.out*cov.p[i,j] + 1.96*deltamethod(as.formula(form),c(c0.out,c1.out),VC[(S+1):(S+2),(S+1):(S+2)]))
				}
			}
					  }

N.low <- exp(outvect[1:S]-1.96*sqrt(diag(VC)[1:S]));N.up <- exp(outvect[1:S]+1.96*sqrt(diag(VC)[1:S])); 



postscript("C://Users/ed234/Dropbox/so_mix/figures/p_temp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(seq(min(data_samp$TEMP,na.rm=TRUE),max(data_samp$TEMP,na.rm=TRUE),.01)*wtmean,expo(fit3$c0.est + fit3$c1.est*seq(min(data_samp$TEMP,na.rm=TRUE),max(data_samp$TEMP,na.rm=TRUE),.01)),xlab="Temperature (Degrees C)",ylab=expression(p),type="l",cex.axis=1.1,cex.lab=1.4)
dev.off()



# standard errors for table
unt_ses <- sqrt(diag(solve(fit1$Hessian)))
unt_ses[51:59]


postscript("P:/Documents/MixtCounts/Paper/Figures/wlat_temp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit3$w.est[,1][order(cov.w)],type="l",ylim=c(0,1),xlab="Northing (km)",ylab=expression(w),cex.axis=1.1,cex.lab=1.4)
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit3$w1.low[order(cov.w)],lty=2,type="l")
lines(unique(data_samp[,c("SITE","NORTH")])$NORTH[order(cov.w)]/1000,fit3$w1.up[order(cov.w)],lty=2,type="l")
dev.off()





postscript("P:/Documents/MixtCounts/Paper/Figures/phitime.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(fit3$phi.est[1,,1],type="l",xlab="Week",ylab=expression(phi),ylim=c(0,1),cex.axis=1.1,cex.lab=1.4)
lines(fit3$phi.low[1,,1],type="l",lty=2)
lines(fit3$phi.up[1,,1],type="l",lty=2)
dev.off()

postscript("P:/Documents/MixtCounts/Paper/Figures/logitphitime.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(logit(fit3$phi.est[1,,1]),type="l",xlab="Week",ylab=expression(logit (phi)),ylim=c(-10,3),cex.axis=1.1,cex.lab=1.4)
lines(logit(fit3$phi.low[1,,1]),type="l",lty=2)
lines(logit(fit3$phi.up[1,,1]),type="l",lty=2)
dev.off()






	
	
	

FittedVal <- FittedVal.f(fit3)
# Residual deviance
ResDev <- ResDev.f(FittedVal)
# DOF
DOF <- length(c(counts[counts > -1]))-fit3$npar
# Residuals

ResDev/DOF


# Exclude zero counts
counts2 <- counts; FittedVal2 <- FittedVal
#counts2[counts2==0] <- -1
FittedVal2[counts2==-1] <- NA
counts2[counts2==-1] <- NA


GOF <- sum((counts2[!is.na(counts2)]-FittedVal2[!is.na(FittedVal2)])^2/FittedVal2[!is.na(FittedVal2)])




ResDev/DOF

postscript("P:/Documents/MixtCounts/Paper/Figures/Ncomp.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
library(plotrix)
plotCI((fit3$N.est),ui=fit3$N.up,li=fit3$N.low,ylab="N",xlab="Site",cex.axis=1.1,cex.lab=1.4)
points(fit1$N.est, pch=19)

#load(file=paste("C://Users/",comp,"/Dropbox/so_mix/weeks/latestCI/fit_",i_species,"_2010_s",S,"_6_SO.RData",sep=""))
#plotCI(fit1$N.est,ui=fit1$N.up,li=fit1$N.low,add=TRUE,cex.axis=1.1,cex.lab=1.4)
dev.off()

c(fit3$inter.est-1.96*sqrt(diag(solve(fit3$Hessian)))[57]*sqrt(ResDev/DOF),fit3$inter.est+1.96*sqrt(diag(solve(fit3$Hessian)))[57])
fit3$slope.est
fit3$slope2.est

round(c(fit3$inter.est,sqrt(diag(solve(fit3$Hessian)))[57]*sqrt(ResDev/DOF)),2)
round(c(fit3$slope.est,sqrt(diag(solve(fit3$Hessian)))[58]*sqrt(ResDev/DOF)),2)

postscript("C://Users/ed234/Dropbox/so_mix/figures/figures_temp/count_fitted.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(log(counts2[counts2 > 0 & !is.na(counts2)]),log(FittedVal2[counts2 > 0 & !is.na(counts2)]),xlab="Log(count)",ylab="Log(fitted)",cex.axis=1.1,cex.lab=1.4)
dev.off()


postscript("C://Users/ed234/Dropbox/so_mix/figures/figures_temp/fitted_res.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 6)
par(mar=c(5.1, 4.4, 4.1, 5.1))
plot(FittedVal,StDevResid.f(FittedVal),xlab="Fitted",ylab="Residuals",cex.lab=1.4,cex.axis=1.1)
dev.off()


















	
	

	
	

