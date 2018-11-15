########################################################
###  STATISTICAL MODELING OF JOINT EXTREMES - R CODE ###
########################################################

### LOADING REQUIRED PACKAGES 
library(ismev)
library(extRemes)
library(evd)
library(extremis)
library(copula)
library(SpatialADAI)

### SLIDE 1- Fox river data
data(fox)
head(fox)
plot(fox,pch=19,main="maximum annual flood discharges (1000 cubed feet per second)")


### SLIDE 2 - Wave and Surge heights
data(wavesurge)
plot(wavesurge,pch=19,cex.lab=1.4)


### SLIDE 7 - Transforming to common scale
Uni_transform_GPD <- function(x,threshold=quantile(x,0.95)){
  par <- fevd(x,threshold = threshold, type="GP")$results$par
  return(ifelse(x>rep(threshold,length(x)),1-length(x[x>threshold])/length(x)*(1+(par[2]/par[1])*(x-threshold))^(-1/par[2]),rank(x)/(length(x)+1)))
}

Frechet_transform <- function(x,threshold=quantile(x,0.95)){
  -1/log(Uni_transform_GPD(x,threshold))
}

GEV2Frechet <- function(x,u,sigma,xi){
  return((1+(xi/sigma)*(x-u))^(1/xi))
}

## SLIDE 8 - Normalized plots
wavesurge_trans <- data.frame(WaveHeight = Frechet_transform(wavesurge[,1]),WaveSurge = Frechet_transform(wavesurge[,2]))
plot(wavesurge_trans,log='xy',pch=19,cex=0.6)

#estimating marginals for Fox data
mar1 <- fevd(fox$berlin,type="GEV")$results$par
mar2 <- fevd(fox$wright,type="GEV")$results$par
plot(GEV2Frechet(fox[,1],mar1[1],mar1[2],mar1[3]),GEV2Frechet(fox[,2],mar2[1],mar2[2],mar2[3]),pch=19,xlab='berlin',ylab='wright')

## SLIDE 10 - Fitting bivariate block-maxima
fit_log <- fbvevd(fox,model="log")
fit_log

## SLIDE 11 - Diagnostics
plot(fbvevd(fox,model="log"))

## SLIDE 12 - Model choice
fit_bilog <- fbvevd(fox,model="bilog")
fit_ct <- fbvevd(fox,model="ct")
dev <- 2*(fit_bilog$deviance-fit_log$deviance)
1-pchisq(dev,1)

## SLIDE 13 - Pickands dependence function
abvnonpar(data=fox,plot=TRUE,lwd=1.5)
abvnonpar(data=fox,add=TRUE,method="pickands",lwd=1.5,lty=2)
abvevd(dep=fit_log$estimate[7],add=TRUE,col="red",lwd=1.5)
abvevd(alpha=fit_bilog$estimate[7],beta=fit_bilog$estimate[8],model="bilog",add=TRUE,col="blue",lwd=1.5)
abvevd(alpha=fit_ct$estimate[7],beta=fit_ct$estimate[8],model="ct",add=TRUE,col="green",lwd=1.5)
legend("bottomleft",col=c("black","black","red","blue","green"),legend=c("emp:cfg","emp:pick","log","bilog","ct"),lty=c(1,2,1,1,1))

## SLIDE 15 - Fitting extreme value copulae
tev <- fitCopula(tevCopula(),pobs(fox))
gumb <- fitCopula(gumbelCopula(),pobs(fox))
hr <- fitCopula(huslerReissCopula(),pobs(fox))
summary(tev)

## SLIDE 16 - Test of goodness of fit
gofEVCopula(gumbelCopula(),fox)
gofEVCopula(tevCopula(),fox,N=100)
gofEVCopula(huslerReissCopula(),fox)

## SLIDE 17 - Pickands dependence function
abvnonpar(data=fox,plot=TRUE,lwd=1.5)
abvnonpar(data=fox,add=TRUE,method="pickands",lwd=1.5,lty=2)
lines(seq(0,1,0.001),A(gumbelCopula(gumb@copula@parameters),seq(0,1,0.001)),col="red",lwd=1.5)
lines(seq(0,1,0.001),A(tevCopula(tev@copula@parameters[1],tev@copula@parameters[2]),seq(0,1,0.001)),col="blue",lwd=1.5)
lines(seq(0,1,0.001),A(huslerReissCopula(hr@copula@parameters[1]),seq(0,1,0.001)),col="green",lwd=1.5)
legend("bottomleft",col=c("black","black","red","blue","green"),legend=c("emp:cfg","emp:pick","gumbel","t-EV","hr"),lty=c(1,2,1,1,1))

## SLIDE 18 - Bivariate exceedances
plot(wavesurge_trans,pch=19,cex.lab=1.4,log='xy')
abline(h=quantile(wavesurge_trans$WaveSurge,0.95),lwd=2)
abline(v=quantile(wavesurge_trans$WaveHeight,0.95),lwd=2)

## SLIDE 20 - Fitting bivariate exceedances
fit1 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)))
fit1

## SLIDE 21 - Model choice
fit2 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="bilog")
fit3 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="ct")

fit1poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),likelihood="poisson")
fit2poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="bilog",likelihood="poisson")
fit3poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="ct",likelihood="poisson")

## SLIDE 22 - Diagnostics
plot(fit1)

## SLIDE 23 - Pickands dependence
abvnonpar(data=wavesurge, method="pot",plot=TRUE,main="Censored")
abvevd(dep=fit1$estimate[5],add=TRUE,col="red")
abvevd(alpha=fit2$estimate[5],beta= fit2$estimate[6],add=TRUE,col="blue",model="bilog")
abvevd(alpha=fit3$estimate[5],beta= fit3$estimate[6],add=TRUE,col="green",model="ct")
legend("bottomleft",col=c("black","red","blue","green"),legend=c("emp","log","bilog","ct"),lty=rep(1,4))

abvnonpar(data=wavesurge, method="pot",plot=TRUE,main="Poisson")
abvevd(dep=fit1poi$estimate[5],add=TRUE,col="red")
abvevd(alpha=fit2poi$estimate[5],beta= fit2poi$estimate[6],add=TRUE,col="blue",model="bilog")
abvevd(alpha=fit3poi$estimate[5],beta= fit3poi$estimate[6],add=TRUE,col="green",model="ct")
legend("bottomleft",col=c("black","red","blue","green"),legend=c("emp","log","bilog","ct"),lty=rep(1,4))

## SLIDE 24 - Spectral density
radius <- (apply(wavesurge_trans,1,sum))
w<- (wavesurge_trans[,1]/radius)[radius >quantile(radius,0.95)]

hist(w,freq=FALSE,main="Censored")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),dep=fit1$estimate[5],half=TRUE),col="red")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit2$estimate[5],beta=fit2$estimate[6],half=TRUE,model="bilog"),col="blue")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit3$estimate[5],beta=fit3$estimate[6],half=TRUE,model="ct"),col="green")

hist(w,freq=FALSE,main="Poisson")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),dep=fit1poi$estimate[5],half=TRUE),col="red")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit2poi$estimate[5],beta=fit2poi$estimate[6],half=TRUE,model="bilog"),col="blue")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit3poi$estimate[5],beta=fit3poi$estimate[6],half=TRUE,model="ct"),col="green")

## SLIDE 25 - River flows
river <- readRDS("river.RDS")
river_trans <- data.frame(Fajardo = Frechet_transform(river[,1]),Espiritu.Santu = Frechet_transform(river[,2]))
plot(river_trans,log='xy',pch=19,cex=0.6)
abline(h=quantile(river_trans$Espiritu.Santu,0.85),lwd=2)
abline(v=quantile(river_trans$Fajardo,0.85),lwd=2)

## SLIDE 26 - Spectral density
fit1 <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)))
fit2 <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),model="bilog")
fit3 <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),model="ct")

fit1poi <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),likelihood = "poisson")
fit2poi <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),model="bilog",likelihood = "poisson")
fit3poi <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),model="ct",likelihood = "poisson")

radius <- (apply(river_trans,1,sum))
w<- (river_trans[,1]/radius)[radius >quantile(radius,0.85)]
hist(w,freq=FALSE,main="Censored",breaks=10)
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),dep=fit1$estimate[5],half=TRUE),col="red")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit2$estimate[5],beta=fit2$estimate[6],half=TRUE,model="bilog"),col="blue")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit3$estimate[5],beta=fit3$estimate[6],half=TRUE,model="ct"),col="green")

hist(w,freq=FALSE,main="Poisson")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),dep=fit1poi$estimate[5],half=TRUE),col="red")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit2poi$estimate[5],beta=fit2poi$estimate[6],half=TRUE,model="bilog"),col="blue")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit3poi$estimate[5],beta=fit3poi$estimate[6],half=TRUE,model="ct"),col="green")


## SLIDE 27 - A different threshold strategy
radius <- (apply(wavesurge_trans,1,sum))
w<- (wavesurge_trans[,1]/radius)[radius >quantile(radius,0.95)]
plot(wavesurge_trans,log='xy',pch=19,cex=0.6)
x<-seq(0.01,quantile(radius,0.95)-0.001,length.out = 100)
y<- quantile(radius,0.95)-x
lines(x,y,lwd=2,col="red")

## SLIDE 29 - Empirical angular measure
radius <- (apply(wavesurge_trans,1,sum))
w<- (wavesurge_trans[,1]/radius)
hist(w[radius >quantile(radius,0.95)],freq=F)

## SLIDE 30 - Non-parametric angular estimation
fit1 <- angdensity(wavesurge_trans, nu = 5, raw = FALSE)
fit2 <- angdensity(wavesurge_trans, nu = 10, raw = FALSE)
fit3 <- angdensity(wavesurge_trans, nu = 50, raw = FALSE)
fit4 <- angdensity(wavesurge_trans, nu = 100, raw = FALSE)
hist(w[radius >quantile(radius,0.95)],freq=F,main="Wavesurge")
rug(w[radius >quantile(radius,0.95)])
lines(fit1$grid,fit1$h,lwd=1.5)
lines(fit2$grid,fit2$h,lwd=1.5,col="red")
lines(fit3$grid,fit3$h,lwd=1.5,col="blue")
lines(fit4$grid,fit4$h,lwd=1.5,col="green")

radius <- (apply(river_trans,1,sum))
w<- (river_trans[,1]/radius)[radius >quantile(radius,0.95)]
fit1 <- angdensity(river_trans, nu = 5, raw = FALSE)
fit2 <- angdensity(river_trans, nu = 10, raw = FALSE)
fit3 <- angdensity(river_trans, nu = 50, raw = FALSE)
fit4 <- angdensity(river_trans, nu = 100, raw = FALSE)
hist(w,freq=F,main="River")
rug(w)
lines(fit1$grid,fit1$h,lwd=1.5)
lines(fit2$grid,fit2$h,lwd=1.5,col="red")
lines(fit3$grid,fit3$h,lwd=1.5,col="blue")
lines(fit4$grid,fit4$h,lwd=1.5,col="green")


## SLIDE 31 - Non-parametric angulard cdf
radius <- (apply(wavesurge_trans,1,sum))
w<- (wavesurge_trans[,1]/radius)[radius >quantile(radius,0.95)]
fit <- angcdf(wavesurge_trans,  raw = FALSE)
fit1 <- angscdf(wavesurge_trans, nu = 5, raw = FALSE)
fit2 <- angscdf(wavesurge_trans, nu = 10, raw = FALSE)
fit3 <- angscdf(wavesurge_trans, nu = 50, raw = FALSE)
fit4 <- angscdf(wavesurge_trans, nu = 100, raw = FALSE)

plot(fit,lwd=2)
lines(fit1$grid,fit1$H,lwd=1.5)
lines(fit2$grid,fit2$H,lwd=1.5,col="red")
lines(fit3$grid,fit3$H,lwd=1.5,col="blue")
lines(fit4$grid,fit4$H,lwd=1.5,col="green")

radius <- (apply(river_trans,1,sum))
w<- (river_trans[,1]/radius)[radius >quantile(radius,0.95)]
fit <- angcdf(river_trans,  raw = FALSE)
fit1 <- angscdf(river_trans, nu = 5, raw = FALSE)
fit2 <- angscdf(river_trans, nu = 10, raw = FALSE)
fit3 <- angscdf(river_trans, nu = 50, raw = FALSE)
fit4 <- angscdf(river_trans, nu = 100, raw = FALSE)

plot(fit,lwd=2)
lines(fit1$grid,fit1$H,lwd=1.5)
lines(fit2$grid,fit2$H,lwd=1.5,col="red")
lines(fit3$grid,fit3$H,lwd=1.5,col="blue")
lines(fit4$grid,fit4$H,lwd=1.5,col="green")

# SLIDE 33 - Threshold diagnostics
bvtcplot(river,type="l", main = "River")
bvtcplot(river,type="l", main = "Wavesurge")

# SLIDE 35 - Chiplots
chiplot(river,nq=400,which=1,spcases=TRUE,ylim1=c(0,1),main1 ="Chi Plot River",xlim=c(0.5,1))
chiplot(wavesurge,nq=400,which=1,spcases=TRUE,ylim1=c(0,1),main1 ="Chi Plot Wavesurge",xlim=c(0.5,1))

# SLIDE 36 - Estimating xhi
2*(1-abvnonpar(0.5,wavesurge))
2*(1-abvnonpar(0.5,river))

fit1 <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)))
fit1poi <- fbvpot(river,threshold = c(quantile(river$Fajardo,0.85),quantile(river$Fajardo,0.85)),likelihood = "poisson")
fit1$dep.summary
fit1poi$dep.summary

fit1 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)))
fit1poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),likelihood="poisson")
fit1$dep.summary
fit1poi$dep.summary

# SLIDE 37 - Chibarplots
chiplot(river,nq=400,which=2,spcases=TRUE,ylim2=c(0,1),main2 ="Chibar Plot River",xlim=c(0.5,1))
chiplot(wavesurge,nq=400,which=2,spcases=TRUE,ylim2=c(0,1),main2 ="Chibar Plot Wavesurge",xlim=c(0.5,1))


fit1 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)))
fit2 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="bilog")
fit3 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="aneglog")
fit4 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="ct")

plot(fit1)

fit1poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),likelihood="poisson")
fit2poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="bilog",likelihood="poisson")
#fit3 <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="aneglog",likelihood="poisson")
fit4poi <- fbvpot(wavesurge,threshold = c(quantile(wavesurge$wave,0.95),quantile(wavesurge$surge,0.95)),model="ct",likelihood="poisson")


## SLIDE 41- Estimating the delta parameter
wavesurge_uni <- cbind(Uni_transform_GPD(wavesurge[,1]),Uni_transform_GPD(wavesurge[,2]))
wv_ievl <- fit.cop.2dim(wavesurge_uni,tresh = 0.95,init.val = c(0.5,0.6),model="IEVL")
wv_gauss <- fit.cop.2dim(wavesurge_uni,0.95,c(0.5,0.6),model="Gauss")

wv_ievl_par <- fit.cop.2dim(wavesurge_uni,0.95,c(0.5,0.6),model="IEVL",cens.type = "partial")
wv_gauss_par0.9 <- fit.cop.2dim(wavesurge_uni,0.9,c(0.5,0.6),model="Gauss",cens.type="partial")

river_uni <- cbind(Uni_transform_GPD(river[,1]),Uni_transform_GPD(river[,2]))
riv_ievl <- fit.cop.2dim(river_uni,0.95,c(0.5,0.6),model="IEVL")
riv_gauss0.9 <- fit.cop.2dim(river_uni,0.9,c(0.5,0.6),model="Gauss")

riv_ievl_par <- fit.cop.2dim(river_uni,0.95,c(0.4,0.3),model="IEVL",cens.type = "partial")
riv_gauss_par0.9 <- fit.cop.2dim(river_uni,0.9,c(0.5,0.6),model="Gauss",cens.type="partial")

## SLIDE 42 - Data on the uniform scale
plot(wavesurge_uni,xlab="wave",ylab="surge")
plot(river_uni,xlab="Fajardo",ylab="Espiritu Santu")

## SLIDE 43 - Coefficient of asymptotic dependence
chiplot(river_uni,qlim=c(0.85,0.99),which=1,xlim=c(0.85,0.99),lwd=2,ylim1=c(0,1), main1="River")
lines(seq(0.8,0.99,0.01),chiu2(seq(0.8,0.99,0.01),riv_ievl$mle[1],riv_ievl$mle[2],"IEVL"),col="red",lwd=2)
lines(seq(0.8,0.99,0.01),chiu2(seq(0.8,0.99,0.01),riv_gauss0.9$mle[1],riv_gauss0.9$mle[2],"Gauss"),col="blue",lwd=2)
lines(seq(0.85,0.99,0.01),chiu2(seq(0.85,0.99,0.01),riv_gauss_par0.9$mle[1],riv_gauss_par0.9$mle[2],"Gauss"),col="green",lwd=2)
legend("bottomleft",col=c("black","red","blue","green"), lty=c(1,1,1,1),legend=c("emp","IEVL","Gauss(0.9)","Gauss_par(0.9)"))

chiplot(wavesurge_uni,qlim=c(0.85,0.99),which=1,xlim=c(0.85,0.99),lwd=2,ylim1=c(0,1),main1="Wavesurge")
lines(seq(0.8,0.99,0.01),chiu2(seq(0.8,0.99,0.01),wv_ievl$mle[1],wv_ievl$mle[2],"IEVL"),col="red",lwd=2)
lines(seq(0.8,0.99,0.01),chiu2(seq(0.8,0.99,0.01),wv_gauss$mle[1],wv_gauss$mle[2],"Gauss"),col="blue",lwd=2)
legend("topleft",col=c("black","red","blue"), lty=c(1,1,1),legend=c("emp","IEVL","Gauss"))
