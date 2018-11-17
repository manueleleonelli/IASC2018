library(tidyverse)
library(extremis)

## SLIDE 5 - Data info
info <- read.csv("info_dep.csv")
head(info)

## SLIDE 6 - Plot Data
rain <- read.csv("Beatenberg_field.csv")
rain <- rain[,1:2]
names(rain) <- c("Date","Temp")
plot(rain$Temp,type="l",xaxt='n',ylab="Temperature",xlab="Date",main = "Beatenberg open field")
a <- seq(1,nrow(rain),length.out=10)
axis(side =1, at=a,rain$Date[a])

## SLIDE 7 - Computing daily maxima
rain1 <- aggregate(Temp ~ format(as.Date(Date), "%Y-%m-%d"), rain, max)
names(rain1) <- c("Date","Temp")
rain1 <- subset(rain1, format(as.Date(Date), "%m-%d") != "02-29")
rain1 <- subset(rain1, format(as.Date(Date), "%Y") != "1997")
rain1 <- subset(rain1, format(as.Date(Date), "%Y") != "2016")
rain1$Date <- as.Date(rain1$Date)
full.dates <- as.Date(as.Date("1998-01-01"):as.Date("2015-12-31"), origin = "1970-01-01")
rain.new <- data.frame(Date = full.dates)
rain.new <- left_join(rain.new,rain1)
rain.new <- subset(rain.new, format(as.Date(Date), "%m-%d") != "02-29")
plot(rain.new$Temp,type="l",xaxt='n',ylab="Temperature",xlab="Date",main = "Beatenberg open field")
a <- seq(1,nrow(rain.new),length.out=10)
axis(side =1, at=a,rain.new$Date[a])

## SLIDE 8 - Fitting a spline
rain2 <- aggregate(Temp ~ format(as.Date(Date), "%m-%d"), rain.new, mean)
names(rain2) <- c("Date","Temp")
x<- 1:365
pispl <- smooth.spline(x , rain2[,2]) 
res <- cbind(Date=format(seq.Date(as.Date("1998-01-01"),by=1,length=365),"%m-%d"),mean=as.numeric(pispl$y))
plot(pispl,type="l",ylim=c(0,20),xlab= "Day of the year", ylab="Average temperature")
points(rain2[,2])

## SLIDE 9 - Normalizing the mean
rain.new$res1 <- rain.new[,2]-pispl$y
plot(rain.new$res1,type="l", xlab="Date",ylab="Normalized temperatures",xaxt='n')
a <- seq(1,nrow(rain.new),length.out=10)
axis(side =1, at=a,rain.new$Date[a])

## SLIDE 10 - Normalizing the variance
rain3 <- aggregate(Temp ~ format(as.Date(Date), "%m-%d"), rain.new, var)
pispl1 <- smooth.spline(x,rain3[,2])
rain.new$res2 <- rain.new[,3]/sqrt(pispl1$y)
plot(rain.new$res2,type='l',xlab="Date",ylab="Normalized temperatures",xaxt='n')
a <- seq(1,nrow(rain.new),length.out=10)
axis(side =1, at=a,rain.new$Date[a])

## SLIDE 12 - Plotting the data
data <- read.csv("Full_data_pairs.csv")
temp <- na.omit(data[,2:3])
temp1 <- na.omit(data[,6:7])
plot(temp,pch=19,cex=0.5,main="Beatenberg")
plot(temp1,pch=19,cex=0.5,main="Celerina")

## SLIDE 13 - Plotting on the uniform scale
temp <-apply(temp,2, rank)/(nrow(temp)+1)
temp1 <-apply(temp1,2, rank)/(nrow(temp1)+1)
plot(temp,pch=19,cex=0.5,main="Beatenberg")
plot(temp1,pch=19,cex=0.5,main="Celerina")

## SLIDE 14 - Plotting on the Frechet scale
Beatenberg <- data.frame(Date=na.omit(data[,c(1:3)])$Date,-1/log(temp))
Celerina <- data.frame(Date=na.omit(data[,c(1,6:7)])$Date,-1/log(temp1))
plot(Beatenberg[,2:3],log='xy',main="Beatenberg",pch=19,cex=0.5)
plot(Celerina[,2:3],log='xy',main="Celerina",pch=19,cex=0.5)

## SLIDE 15 - Standard modeling
Beat_radius <- apply(Beatenberg[,2:3],1,sum)
Cel_radius <- apply(Celerina[,2:3],1,sum)
plot(Beatenberg[,2:3],log='xy',main="Beatenberg",pch=19,cex=0.5)
x<-seq(0.01,quantile(Beat_radius,0.98)-0.001,length.out = 100)
y<- quantile(Beat_radius,0.98)-x
lines(x,y,lwd=2,col="red")
plot(Celerina[,2:3],log='xy',main="Celerina",pch=19,cex=0.5)
x<-seq(0.01,quantile(Cel_radius,0.98)-0.001,length.out = 100)
y<- quantile(Cel_radius,0.98)-x
lines(x,y,lwd=2,col="red")

## SLIDE 16 - Angular densities
Beat_w <- (Beatenberg[,2]/Beat_radius)[Beat_radius>quantile(Beat_radius,0.98)]
Cel_w <- (Celerina[,2]/Cel_radius)[Cel_radius>quantile(Cel_radius,0.98)]
hist(Beat_w,freq=F,main="Beatenberg",ylim=c(0,2))
Beat_dens<-angdensity(Beatenberg[,2:3],tau=0.98,100,raw=F)
lines(Beat_dens$grid,Beat_dens$h,lwd=2,col='red')
hist(Cel_w,freq=F,main="Celerina",ylim=c(0,2))
Cel_dens<-angdensity(Celerina[,2:3],tau=0.98,100,raw=F)
lines(Cel_dens$grid,Cel_dens$h,lwd=2,col='red')



ext <- read.csv("canopy_dependence.csv")
head(ext)


## IMPLEMENTATION OF THE SERR ESTIMATOR
serr <- function(d,nu,b,grid_w,grid_x){
  N <- nlevels(as.factor(d[,2])) 
  x <- as.numeric(levels(as.factor(d[,2])))
  temp <- split(d[,1],d[,2])
  dens <- matrix(0,length(grid_w),N)
  for(j in 1:N) {
    dens[,j] <- angdens(temp[[j]],nu,grid_w)
  }
  hat_h <- matrix(0,length(grid_w),length(grid_x))
  for(i in 1:length(grid_w)){
    for(j in 1:length(grid_x)){
      hat_h[i,j] <- sum(dnorm((grid_x[j]-x)/b)*dens[i,])/sum(dnorm((grid_x[j]-x)/b))
    }
  }
  return(hat_h)
}

angdens <- function(d,nu,grid_w){
  temp <- rep(0,length(grid_w))
  for (i in 1:length(grid_w)){
    temp[i] <- sum(p(d)*dbeta(grid_w[i],nu*d,nu*(1-d)))
  }
  temp
}

p <- function(d){
  (1/length(d))*(1-(mean(d)-0.5)*var(d)*(d-mean(d)))
}


## ALTITUDE
grid_w=seq(0.01,0.99,length.out=100)
grid_x =seq(min(ext$Altitude),max(ext$Altitude),length.out=100)
fit <- serr(ext[,2:3],10,50,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- serr(ext[,2:3],50,50,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- serr(ext[,2:3],50,100,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- serr(ext[,2:3],10,100,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")


## IMPLEMENTATION OF THE AoAS ESTIMATOR
AoAS <-function(w,x,nu,b,tau,grid_w,grid_x){
  h<- matrix(0,length(grid_w),length(grid_x))
  for(i in 1:length(grid_w)){
    for(j in 1:length(grid_x)){
      pi <- (1/b)*dnorm((grid_x[j]-x)/b)/sum((1/b)*dnorm((grid_x[j]-x)/b))
      theta <- 0.5/sum(pi*w)
      h[i,j] <- sum(pi*dbeta(grid_w[i],nu*w*theta+tau,nu*(1-w*theta)+tau))
    }
  }
  return(h)
}

## SLOPE
grid_w=seq(0.01,0.99,length.out=100)
grid_x =seq(min(ext$Slope),max(ext$Slope),length.out=100)
fit <- AoAS(ext[,2],ext[,6],100,5,1,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- AoAS(ext[,2],ext[,6],20,5,1,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- AoAS(ext[,2],ext[,6],100,20,1,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

fit <- AoAS(ext[,2],ext[,6],100,10,1,grid_w,grid_x)
persp3D(grid_w,grid_x,fit,xlab="w",ylab="x",zlab="dens")

