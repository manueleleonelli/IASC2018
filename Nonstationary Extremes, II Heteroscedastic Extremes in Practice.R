#####################################################
###  NON-STATIONARY EXTREMES IN PRACTICE - R CODE ###
#####################################################

### LOADING REQUIRED PACKAGES 
library(ismev)
library(extRemes)
library(extremis)
library(evmix)
library(ggplot2)

## SLIDE 5 - Plotting Fremantle data
data(fremantle)
head(fremantle)
plot(fremantle[,1:2],xlab = "Year",ylab="Sea Level in ms",pch=19)


## SLIDE 6 - Plotting Sea Level vs Southern Oscillation Index
sp2<-ggplot(fremantle, aes(x=SOI, y=SeaLevel, color=Year)) + geom_point()
sp2

## SLIDE 8 - Fitting EVD distributions
fremantle_fit <- fevd(fremantle$SeaLevel,type="GEV")
summary(fremantle_fit)

## SLIDE 9 - QQplot
plot(fremantle_fit,type = "qq")

## SLIDE 10 - Density Plot
plot(fremantle_fit,type = "density")

## SLIDE 11 - Return Plot
plot(fremantle_fit,type = "rl")

## SLIDE 12 - Confidence Intervals
ci(fremantle_fit,type="parameter")

## SLIDE 13 - Return levels
return.level(fremantle_fit,return.period=c(2,10,20,50,100,500),do.ci=TRUE)

## SLIDE 14 - Bayesian GEV fitting
fremantle_fit_Bayes <- fevd(fremantle$SeaLevel,type="GEV",method="Bayesian")
summary(fremantle_fit_Bayes)

## SLIDE 15 - Bayesian return plots
plot(fremantle_fit_Bayes,"rl",main="")

## SLIDE 16 - Convergence checking
plot(fremantle_fit_Bayes,type="trace")

## SLIDE 17 - Estimated location
plot(fremantle[,1:2],xlab = "Year",ylab="Sea Level in ms",pch=19)
abline(h=fremantle_fit$results$par[1])

## SLIDE 18 - Non-stationary fitting via time covariate
fremantle$Time <- fremantle$Year-min(fremantle$Year)
fremantle_fit_time <- fevd(SeaLevel,data=fremantle, 
                           location.fun = ~Time)
summary(fremantle_fit_time)

## SLIDE 19 - Is there a time effect?
ci(fremantle_fit_time,type="parameter")

## SLIDE 20 - Time varying return periods
plot(fremantle_fit_time,"rl",rperiods=c(5,10,50))

## SLIDE 21 - Time varying location parameter
plot(fremantle[,4], fremantle[,2],xlab = "Year",ylab="Sea Level in ms",pch=19,xaxt='n')
axis(side =1, at=fremantle$Time[seq(1,nrow(fremantle),by=4)],fremantle$Year[seq(1,nrow(fremantle),by=4)])
abline(a=fremantle_fit_time$results$par[1],b=fremantle_fit_time$results$par[2])

## SLIDE 22 - Model selection via hypothesis testing
dev <- 2*(-fremantle_fit_time$results$value+fremantle_fit$results$value)
dev > qchisq(0.95,1)
1-pchisq(dev,1)

## SLIDE 23 - Quadratic time effect
fremantle_fit_time2 <- fevd(SeaLevel,data=fremantle, location.fun = ~Time + I(Time^2))
plot(fremantle[,4], fremantle[,2],xlab = "Year",ylab="Sea Level in ms",pch=19,xaxt='n')
axis(side =1, at=fremantle$Time[seq(1,nrow(fremantle),by=4)],fremantle$Year[seq(1,nrow(fremantle),by=4)])
lines(fremantle$Time,fremantle_fit_time2$results$par[1]+fremantle_fit_time2$results$par[2]*fremantle$Time+fremantle_fit_time2$results$par[3]*fremantle$Time^2)

dev <- 2*(-fremantle_fit_time2$results$value+fremantle_fit_time$results$value)
1- pchisq(dev,1)

## SLIDE 24 - Adding SOI to the model
fremantle_fit_SOI <- fevd(SeaLevel,data=fremantle, location.fun = ~SOI+Time)
ci(fremantle_fit_SOI,type="parameter")

## SLIDE 25 - Return levels depending on SOI and time
plot(fremantle_fit_SOI,type="rl",rperiod=c(5,10,50))
dev <- 2*(-fremantle_fit_SOI$results$value+fremantle_fit$results$value)
1- pchisq(dev,1)

## SLIDE 26 - Time/SOI interaction
fremantle_fit_SOItime <- fevd(SeaLevel,data=fremantle, location.fun = ~SOI*Time)
sp2<-ggplot(fremantle, aes(x=SOI, y=SeaLevel, color=Year)) + geom_point()
sp2

dev <- 2*(-fremantle_fit_SOItime$results$value+fremantle_fit_SOI$results$value)
1-pchisq(dev,1)

## SLIDE 27 - Fort Data
data(Fort)
head(Fort)

## SLIDE 28 - Extracting maxima
bmFort <- blockmaxxer(Fort, blocks = Fort$year, which="Prec")
plot(Fort$year, Fort$Prec, xlab = "Year", ylab = "Precipitation (inches)",
     cex = 1.25, cex.lab = 1.25,
     col = "darkblue", bg = "lightblue", pch = 21)
points(bmFort$year, bmFort$Prec, col="darkred", cex=1.5)

## SLIDE 29 - Venice sea levels
data(venice)
head(venice)

## SLIDE 30 - Venice sea level data
plot(venice[,1:2],ylim=c(min(na.omit(venice[,-1])),max(na.omit(venice[,-1]))),ylab="SeaLevel")
for(i in 3:ncol(venice)) points(venice$Year,venice[,i])

## SLIDE 31 - Fitting r-largest observations
venice <- venice[-5,] #omit NA
venice$time <- seq(1:nrow(venice)) #add time variable
venice_fit_time <-fevd(r1,type="GEV", data = venice,location.fun = ~time)
venice_fit_time_5 <- rlarg.fit(venice[,2:6], ydat=venice,mul=12)
venice_fit_time_10 <- rlarg.fit(venice[,2:11],ydat=venice, mul=12)


## SLIDE 32 - Sea level location
plot(venice[,12],venice[,2],ylim=c(min(na.omit(venice[,-c(1,12)])),max(na.omit(venice[,-c(1,12)]))))
for(i in 3:(ncol(venice)-2)) points(venice$time,venice[,i])
abline(a=venice_fit_time_5$mle[1],b=venice_fit_time_5$mle[2],col="blue")
abline(a=venice_fit_time_10$mle[1],b=venice_fit_time_10$mle[2],col="red")
abline(a=venice_fit_time$results$par[1],b=venice_fit_time$results$par[2])
legend("topleft",col=c("black","blue","red"),lty=c(1,1,1),legend=c("r=1","r=5","r=10"))

## SLIDE 33 - Dow Jones Index
data(dowjones)
plot(dowjones,type="l")

## SLIDE 34 - Negative log-returns
dowjones <- data.frame(dowjones[-1,],returns=-100*diff(log(dowjones$Index)))
plot(dowjones$Date,dowjones$returns,type="l")

## SLIDE 35 - The returns' distribution
hist(dowjones$returns)

## SLIDE 36 - MRL plot
mrl.plot(dowjones$returns)

## SLIDE 37 - A nicer MRL plot
mrlplot(dowjones$returns,try.thresh=c(quantile(dowjones$returns,0.9),2),ylim = c(-2,4))

## SLIDE 38 - Shape stability plot
tshapeplot(dowjones$returns,ylim=c(-2,2),try.thresh = 2)

## SLIDE 39 - Scale stability plot
tscaleplot(dowjones$returns,ylim=c(-8,4),try.thresh = 2)

## SLIDE 40 - Hill plot
hillplot(dowjones$returns,orderl=c(1,150),try.thresh= 2,xlab="Order Statistic")

## SLIDE 41 - Fitting the GPD
dow_fit <- fevd(dowjones$returns,type="GP",threshold=2)
summary(dow_fit)

## SLIDE 42- Histograms
dow_fit2 <- fevd(dowjones$returns,type="GP",threshold=quantile(dowjones$returns,0.9))
plot(dow_fit,"hist")
plot(dow_fit2,"hist")

## SLIDE 43 - Posterior summaries
ci(dow_fit,type="parameter")
ci(dow_fit2,type="parameter")
return.level(dow_fit,do.ci=TRUE)
return.level(dow_fit2,do.ci=TRUE)

## SLIDE 44 - Quantile plots
plot(dow_fit,"qq")
plot(dow_fit2,"qq")

## SLIDE 45 - Are exceedances independent?
dowjones$thresh <- dowjones$returns>2
plot(dowjones$Date,dowjones$returns, xlab="Year",ylab="returns",ylim=c(0,max(dowjones$returns)))
points(subset(dowjones,dowjones$thres==FALSE)$Date,subset(dowjones,dowjones$thres==FALSE)$returns,col="white")
abline(h=2,lwd=2)

dowjones$thresh <- dowjones$returns>quantile(dowjones$returns,0.9)
plot(dowjones$Date,dowjones$returns, xlab="Year",ylab="returns",ylim=c(0,max(dowjones$returns)))
points(subset(dowjones,dowjones$thres==FALSE)$Date,subset(dowjones,dowjones$thres==FALSE)$returns,col="white")
abline(h=quantile(dowjones$returns,0.9),lwd=2)

## SLIDE 46 - Extremal dependence
i <- 1:50
a <- b <- numeric(length(i))
for (j in 1:length(i)){
  a[j] <- extremalindex(dowjones$returns,threshold=2,method="runs",r=j)[1]
  b[j] <- extremalindex(dowjones$returns,threshold=quantile(dowjones$returns,0.9),method="runs",r=j)[1]
}
plot(i,a,xlab="run",ylab="Extremal Index",type="l",lwd=2,ylim=c(0,1))
lines(i,b,lwd=2,col='red')
abline(h=extremalindex(dowjones$returns,threshold=2,method="interval"),lty=2,lwd=2)
abline(h=extremalindex(dowjones$returns,threshold=quantile(dowjones$returns,0.9),method="interval"),lty=2,lwd=2,col='red')
legend("topright",lty=c(1,2,1,2),col=c("black","black","red","red"),legend=c("run/u=2","interval/u=2","run/u=q(0.9)","interval/u=q(0.9)"))


## SLIDE 47 - Declustering
dow_decluster1<-decluster(dowjones$returns,threshold=2,method="runs",r=7,clusterfun = "max")
dow_decluster2<-decluster(dowjones$returns,threshold=quantile(dowjones$returns,0.9),method="interval")
dow_decluster3<-decluster(dowjones$returns,threshold=2,method="interval")
dow_decluster4<-decluster(dowjones$returns,threshold=quantile(dowjones$returns,0.9),method="runs",r=7)

par(mfrow=c(2,2))
plot(dow_decluster3,ylim=c(1,7), main="")
plot(dow_decluster1,ylim=c(1,7),main="")
plot(dow_decluster2,ylim=c(1,7),main="")
plot(dow_decluster4,ylim=c(1,7),main="")

## SLIDE 48 - Results after declustering
fevd(dow_decluster1,threshold=2,type="GP")
fevd(dow_decluster2,threshold=quantile(dowjones$returns,0.9),type="GP")
fevd(dow_decluster3,threshold=2,type="GP")
fevd(dow_decluster4,threshold=quantile(dowjones$returns,0.9),type="GP")


## SLIDE 50 - Results for dow-jones series
# definition of Hill estimator
hill <- function(x,k){
  x <- sort(x)
  (1/k)*sum(log(x[(length(x)-k+1):length(x)]))-log(x[length(x)-k])
}

# Test of Einmhal et al. for constant xi
# m = 4, u =2
m<-4
k <- nrow(subset(dowjones,dowjones$returns>2))
subk<- floor(k/m)+1
gam_par <- c(hill(dowjones$returns[1:325],subk), hill(dowjones$returns[326:650],subk), hill(dowjones$returns[651:925],subk), hill(dowjones$returns[926:nrow(dowjones)],subk))
gam_tot <- hill(dowjones$returns,k)
test<- (1/m)*sum(gam_par/gam_tot-1)^2
1-pchisq(k*test,m-1)

# m = 4 , u = q(0.9)
m<-4
k <- round(nrow(dowjones)*0.1)
subk<- floor(k/m)+1
gam_par <- c(hill(dowjones$returns[1:325],subk), hill(dowjones$returns[326:650],subk), hill(dowjones$returns[651:925],subk), hill(dowjones$returns[926:nrow(dowjones)],subk))
gam_tot <- hill(dowjones$returns,k)
test<- (1/m)*sum(gam_par/gam_tot-1)^2
1-pchisq(k*test,m-1)

# m = 3, u = 2
m <- 3
k <- nrow(subset(dowjones,dowjones$returns>2))
subk<- floor(k/m)+1
gam_par <-  c(hill(dowjones$returns[1:430],subk), hill(dowjones$returns[431:860],subk),hill(dowjones$returns[861:nrow(dowjones)],subk))
gam_tot <- hill(dowjones$returns,k)
test<- (1/m)*sum(gam_par/gam_tot-1)^2
1-pchisq(k*test,m-1)

# m = 3, u = q(0.9)
m <- 3
k <- round(nrow(dowjones)*0.1)
subk<- floor(k/m)+1
gam_par <-  c(hill(dowjones$returns[1:430],subk), hill(dowjones$returns[431:860],subk),hill(dowjones$returns[861:nrow(dowjones)],subk))
gam_tot <- hill(dowjones$returns,k)
test<- (1/m)*sum(gam_par/gam_tot-1)^2
1-pchisq(k*test,m-1)

## SLIDE 51 - Scedasis function
plot(cdensity(dowjones[,-2],2))
plot(cdensity(dowjones[,-2],quantile(dowjones$returns,0.9)))

# SLIDE 52 - Scedasis cdf
plot(cdf(dowjones[,-2],threshold=2))
plot(cdf(dowjones[,-2],threshold=quantile(dowjones$returns,0.9)))

