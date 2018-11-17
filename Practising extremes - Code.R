## Loading the packages ##
library(extRemes)
library(ismev)
library(evmix)
library(evd)
library(copula)
library(extremis)



##### TASK 1 #####
# Loading the data and plotting it
data(ftcanmax)
plot(ftcanmax) #Some extremes, especially in recent years.

# Fitting a GEV
fit <- fevd(ftcanmax[,2],type="GEV")
fit #positive shape estimate but fairly large standard error
plot(fit) #The fit looks reasonable
ci(fit,type="parameter") #As expected the 95% confidence interval includes negative values
return.level(fit,do.ci=T)

#Let's try and see if there is significant time variation of location
ftcanmax$Time <- ftcanmax$Year-min(ftcanmax$Year)
fit_time <- fevd(ftcanmax[,2],data=ftcanmax,location.fun = ~Time)
fit_time #Does not seem to be significant the extra parameter
plot(fit_time)
plot(fit_time,type="rl",rperiods = c(2,5,10,50)) #Return levels are somehow constant

dev <- 2*(-fit_time$results$value + fit$results$value) #Computation of the deviance
1-pchisq(dev,1) #deviance test does not reject the simpler model

##### TASK 2 #####
data(euroex)
plot(euroex,type="l") #Closing prices are not the object of interest
euroex <- -100*diff(log(euroex)) #computing the negative log returns
plot(euroex,type="l")
hist(euroex) #the tail does not seem to be particularly heavy

#Let's try and choose the threshold
mrlplot(euroex)
 
# I picked two possible thresholds
mrlplot(euroex,try.thresh = c(0.5,0.8))
tshapeplot(euroex,try.thresh = c(0.5,0.8))
tscaleplot(euroex,try.thresh = c(0.5,0.8),ylim=c(-2,3))
#Because of the very negative shape, the software gives some warnings
# We cannot use the hillestimator with negative shape

#Fitting the model
fit <- fevd(euroex,type="GP",threshold=0.5)
fit #negative shape estimate
plot(fit) #the fit looks good
return.level(fit) #returns levels do not increase substantially

#Trying a different threshold
#Fitting the model
fit <- fevd(euroex,type="GP",threshold=0.8)
fit #negative shape estimate
plot(fit) #the fit looks good
return.level(fit) #returns levels do not increase substantially
#the results look fairly similar for both thresholds

##### TASK 3 #####
data(wind)
plot(wind[,2:3])
fit <- fbvevd(wind[,2:3],model="log")
fit # high estimate of the dependence parameter - perhaps independence
plot(fit) # fit looks acceptable

#fitting other parametric model
fit_bilog <- fbvevd(wind[,2:3],model="bilog")
fit_bilog
fit_ct <- fbvevd(wind[,2:3],model="ct")
fit_ct

# bilogistic seems to be the best fit
plot(fit_bilog)

#Creating the Pickands dependence approach
abvnonpar(data=wind[,2:3],plot=T)
abvnonpar(data=wind[,2:3],add=T, method="pickands")
abvevd(dep=fit$estimate[7],model="log",add=T,col="blue")
abvevd(alpha=fit_bilog$estimate[7],beta = fit_bilog$estimate[8],model="bilog",add=T,col="red")
abvevd(alpha=fit_ct$estimate[7],beta = fit_ct$estimate[8],model="ct",add=T,col="green")
# actually the Coles-Tawn model seems to provide a better fit.

#Fitting extreme value copulae
tev <- fitCopula(tevCopula(),pobs(wind[,2:3]))
gumb <- fitCopula(gumbelCopula(),pobs(wind[,2:3]))
hr <- fitCopula(huslerReissCopula(),pobs(wind[,2:3]))
summary(tev)
summary(gumb)
summary(hr)
# Pickands dependence function, again the fit looks identical
abvnonpar(data=wind[,2:3],plot=TRUE,lwd=1.5)
abvnonpar(data=wind[,2:3],add=TRUE,method="pickands",lwd=1.5,lty=2)
lines(seq(0,1,0.001),A(gumbelCopula(gumb@copula@parameters),seq(0,1,0.001)),col="red",lwd=1.5)
lines(seq(0,1,0.001),A(tevCopula(tev@copula@parameters[1],tev@copula@parameters[2]),seq(0,1,0.001)),col="blue",lwd=1.5)
lines(seq(0,1,0.001),A(huslerReissCopula(hr@copula@parameters[1]),seq(0,1,0.001)),col="green",lwd=1.5)

##### TASK 4 #####
data(exchange)
exchange <- data.frame(USD.GBP=-100*diff(log(exchange[,1])),CAD.GBP=-100*diff(log(exchange[,2])))
#Plotting the data
plot(exchange,pch=19,cex=0.5)
abline(h=quantile(exchange[,1],0.9),lwd=2)
abline(v=quantile(exchange[,2],0.9),lwd=2)

# Fitting a logistic model
fit1 <- fbvpot(exchange,threshold = c(quantile(exchange[,1],0.9),quantile(exchange[,2],0.9)))
fit1
plot(fit1)

# Fitting other models
fit2 <- fbvpot(exchange,threshold = c(quantile(exchange[,1],0.9),quantile(exchange[,2],0.9)),model="bilog")
fit3 <- fbvpot(exchange,threshold = c(quantile(exchange[,1],0.9),quantile(exchange[,2],0.9)),model="ct")

fit2
fit3
# AIC strongly suggests the bilogistic model

# Producing histogram of the spectral density
# Transforming data into common margins
Uni_transform_GPD <- function(x,threshold=quantile(x,0.95)){
  par <- fevd(x,threshold = threshold, type="GP")$results$par
  return(ifelse(x>rep(threshold,length(x)),1-length(x[x>threshold])/length(x)*(1+(par[2]/par[1])*(x-threshold))^(-1/par[2]),rank(x)/(length(x)+1)))
}

Frechet_transform <- function(x,threshold=quantile(x,0.95)){
  -1/log(Uni_transform_GPD(x,threshold))
}

exchange.new <- data.frame(USD.GBP=Frechet_transform(exchange[,1],0.9),CAD.GBP=Frechet_transform(exchange[,2],0.9))

#Creating the spectral sample
radius <- (apply(exchange.new,1,sum))
w<- (exchange.new[,1]/radius)[radius >quantile(radius,0.9)]

# Plotting the sample and the estimated angular densities
hist(w,freq=F)
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),dep=fit1$estimate[5],half=TRUE),col="red")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit2$estimate[5],beta=fit2$estimate[6],half=TRUE,model="bilog"),col="blue")
lines(seq(0,1,0.01),hbvevd(seq(0,1,0.01),alpha=fit3$estimate[5],beta=fit3$estimate[6],half=TRUE,model="ct"),col="green")
## both the bilogistic and the ct models seem to give a reasonable fit

# Fitting non-parametrically the angular density
ang1 <- angdensity(exchange.new, nu = 5, raw = FALSE)
ang2 <- angdensity(exchange.new, nu = 10, raw = FALSE)
ang3 <- angdensity(exchange.new, nu = 50, raw = FALSE)
ang4 <- angdensity(exchange.new, nu = 100, raw = FALSE)

# Plotting the angular density
hist(w,freq=F)
rug(w)
lines(ang1$grid,ang1$h,lwd=1.5)
lines(ang2$grid,ang2$h,lwd=1.5,col="red")
lines(ang3$grid,ang3$h,lwd=1.5,col="blue")
lines(ang4$grid,ang4$h,lwd=1.5,col="green")

## Assessing extreme dependence
chiplot(exchange,which=1,ylim1=c(0,1))
# the plot does not seem to go to zero, indicating dependence

chiplot(exchange,which=2,ylim2=c(0,1))
# this is confirmed by the chibarplot