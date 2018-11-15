##########################################################
###  STATISTICAL MODELING OF SPATIAL EXTREMES - R CODE ###
##########################################################

library(SpatialExtremes)
library(fields)
library(maps)
library(mapdata)

## SLIDE 3 - Plot the map
par(mar=rep(0,4),ps=16)
image(lon.vec,lat.vec,alt.mat,col=terrain.colors(64),asp=1,bty='n',xlab='',ylab='',axes=F)
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)


## SLIDE 4 - Swiss precipitation dataset
data(rainfall)
data(swissalt)
Year <- seq(1962,2008)

## SLIDE 5 - Spatial dependence
plot(rain[,3],type="l",xaxt='n',xlab="Year",ylab="Rainfall",lwd=1.5)
axis(1,at=seq(1:47),labels=Year)
lines(rain[,36],col="red",lwd=1.5)
plot(rain[,15],type="l",xaxt='n',xlab="Year",ylab="Rainfall",lwd=1.5)
axis(1,at=seq(1:47),labels=Year)
lines(rain[,12],col="red",lwd=1.5)


## SLIDE 8/10 - Exploratory graphs
symbolplot(rain,coord,which="gev",plot.border = swiss,plot.legend=F,scale=2)


## SLIDE 11 - Fitting spatial GEVs
loc.form <- scale.form <- y~lon+lat 
shape.form <- y~1

fit <-fitspatgev(rain,coord,loc.form,scale.form,shape.form)
fit

## SLIDE 12 - Model selection
loc.form <- y~lon*lat
scale.form <- y~ lon+lat
fit1 <- fitspatgev(rain,coord,loc.form,scale.form,shape.form)

loc.form <- shape.form <- scale.form <- y~ 1
fit2 <- fitspatgev(rain,coord,loc.form,scale.form,shape.form)

TIC(fit,fit1,fit2)

## SLIDE 13 - Model selection
anova(fit2,fit)
anova(fit,fit1)

## SLIDE 14 - Return levels
x <- seq(650,800,length=100)
y <- seq(200,300,length=100)
grid <- expand.grid(x,y)
colnames(grid) <- c("lon","lat")
pred <- predict(fit,newdata=grid,ret.per=25)
pred1 <- predict(fit1,newdata=grid,ret.per=25)

image(x,y,matrix(pred$Q25,100),col=tim.colors(64),xlab="Lat",ylab="Lon")
contour(x,y,matrix(pred$Q25,100),col="black",add=TRUE,lwd=2)
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

image(x,y,matrix(pred1$Q25,100),col=tim.colors(64),xlab="Lat",ylab="Lon")
contour(x,y,matrix(pred1$Q25,100),col="black",add=TRUE,lwd=2)
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)


## SLIDES 16/17 - Setting up a latent model
loc.form <- scale.form <- y~lon+lat; shape.form <- y~1

hyper <- list()
hyper$betaMeans <- list(loc=rep(0,3),scale=rep(0,3),shape=0)
hyper$betaIcov <- list(loc= diag(rep(1/1000,3)),scale= diag(rep(1/1000,3)),shape=1/10)

hyper$sills <- list(loc= c(1,12), scale=c(1,1), shape = c(1,0.04))
hyper$ranges <- list(loc = c(5,3), scale= c(5,3), shape = c(5,3))
hyper$smooths <- list(loc = c(1,1), scale= c(1,1), shape = c(1,1))

prop <- list(gev=c(3,0.1,0.3), ranges=c(1,0.8,1.2), smooths=rep(0,3))

start <- list(sills=c(10,10,0.5), ranges=c(20,10,10),smooths=c(1,1,1), beta=list(loc=c(25,0,0),scale=c(33,0,0),shape=0.001))

chain <- latent(rain,coord[,1:2],"powexp",loc.form,scale.form,shape.form,hyper=hyper,prop=prop,start=start,n=1000,burn.in=1000,thin=10)

## SLIDES 18/19 - Model output
chain

## SLIDE 20 - Return levels
map.latent(chain,ret.per=25,plot.contour = T)



## SLIDE 21 - Simulating from the estimated model
post_loc <- apply(chain$chain.loc,2,mean)[1:6]
post_scale <- apply(chain$chain.scale,2,mean)[1:6]
post_shape <- apply(chain$chain.shape,2,mean)[1:3]
n<- 50
n.obs <-1
x <- seq(650,800,length=n)
y <- seq(200,300,length=n)
coord_sim <- cbind(lon= rep(x,n),lat = rep(y,each=n))
gp.loc <- rgp(1, coord_sim, "powexp", sill = post_loc[4], range = post_loc[5], smooth = 1)
gp.scale <- rgp(1, coord_sim, "powexp", sill = post_scale[4], range = post_scale[5], smooth = 1)
gp.shape <- rgp(1, coord_sim, "powexp", sill = post_shape[2], range = post_shape[3], smooth = 1)
locs <- post_loc[1] + post_loc[2] * coord_sim[,"lon"] +post_loc[3]*coord_sim[,"lat"] + gp.loc
scales <- post_scale[1] + post_scale[2] * coord_sim[,"lon"] +post_scale[3]*coord_sim[,"lat"]  + gp.scale
shapes <- post_shape[1]+ gp.shape
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- rgev(n.obs, locs[i], scales[i], shapes[i])
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")


## SLIDE 23 - Fit Gaussian copula
fit_whitmat <- fitcopula(rain,coord[,1:2],copula="gaussian",cov.mod = "whitmat",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_whitmat


## SLIDE 24 - Fit various other copulae
fit_stable <- fitcopula(rain,coord[,1:2],copula="gaussian",cov.mod = "powexp",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_cauchy <- fitcopula(rain,coord[,1:2],copula="gaussian",cov.mod = "cauchy",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_exp <- fitcopula(rain,coord[,1:2],copula="gaussian",cov.mod = "powexp",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form,smooth=1)

fit_tstable <- fitcopula(rain,coord[,1:2],copula="student",cov.mod = "powexp",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_twhitmat <- fitcopula(rain,coord[,1:2],copula="student",cov.mod = "whitmat",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_tcauchy <- fitcopula(rain,coord[,1:2],copula="student",cov.mod = "cauchy",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form)
fit_texp <- fitcopula(rain,coord[,1:2],copula="student",cov.mod = "powexp",loc.form=loc.form,shape.form = shape.form,scale.form=scale.form,smooth=1)

c(AIC(fit_stable),AIC(fit_whitmat),AIC(fit_cauchy),AIC(fit_exp),AIC(fit_tstable),AIC(fit_whitmat),AIC(fit_tcauchy),AIC(fit_texp))


## SLIDE 25 - F-madogram
par(mar=c(4,5,1,1))
fmadogram(rain, coord, fit_whitmat, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, fit_whitmat, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)
fmadogram(rain, coord, fit_texp,n.bin=300, which = "ext",add=T,lty=2,pch=19,cex=0.4,lwd=2)

## SLIDE 26 - 25years return levels
x <- seq(630,780,length=100)
y <- seq(200,300,length=100)
grid <- expand.grid(x,y)
colnames(grid) <- c("lon","lat")
pred_gaus <- predict(fit_whitmat,newdata=grid,ret.per=25)
pred_t <- predict(fit_texp,newdata=grid,ret.per=25)
image(x,y,matrix(pred_t$Q25,100),col=tim.colors(64),xlab="Lon",ylab="Lat")
contour(x,y,matrix(pred_t$Q25,100),col="black",add=TRUE,lwd=2)
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

image(x,y,matrix(pred_gaus$Q25,100),col=tim.colors(64),xlab="Lon",ylab="Lat")
contour(x,y,matrix(pred_t$Q25,100),col="black",add=TRUE,lwd=2)
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)


## SLIDE 27 - Simulations from the estimated models
n<- 75
n.obs <-1
x <- seq(630,780,length=n)
y <- seq(200,300,length=n)
coord_sim <- cbind(lon= rep(x,n),lat = rep(y,each=n))
# simulating the gaussian copula
sim_gaus_cop <- rcopula(1,as.matrix(coord_sim),copula="gaussian",cov.mod="whitmat",smooth=fit_whitmat$param[3],nugget=fit_whitmat$param[1],range=fit_whitmat$param[2])
# define the spatial model for the GEV parameters
param.loc <- fit_whitmat$param[4] + fit_whitmat$param[5] * coord_sim[,1]+ fit_whitmat$param[6]*coord_sim[,2]
param.scale <- fit_whitmat$param[7] + fit_whitmat$param[8] * coord_sim[,1]+ fit_whitmat$param[9]*coord_sim[,2]
param.shape <- fit_whitmat$param[10]
# creating the data
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(sim_gaus_cop[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

# simulating the gaussian copula
sim_t_cop <- rcopula(1,as.matrix(coord_sim),copula="student",cov.mod="powexp",smooth=fit_texp$param[4],nugget=fit_texp$param[2],range=fit_texp$param[3],DoF=fit_texp$param[1])
# define the spatial model for the GEV parameters
param.loc <- fit_texp$param[5] + fit_texp$param[6] * coord_sim[,1]+ fit_texp$param[7]*coord_sim[,2]
param.scale <- fit_texp$param[8] + fit_texp$param[9] * coord_sim[,1]+ fit_texp$param[10]*coord_sim[,2]
param.shape <- fit_texp$param[11]
# creating the data
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(sim_t_cop[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

## SLIDE 28 - Fitting extreme-T copula
ext_stab <- fitmaxstab(rain,coord[,1:2],"tpowexp",loc.form,scale.form,shape.form)
ext_exp <- fitmaxstab(rain,coord[,1:2],"tpowexp",loc.form,scale.form,shape.form,smooth=1)
ext_cauchy <- fitmaxstab(rain,coord[,1:2],"tcauchy",loc.form,scale.form,shape.form)
ext_matern <- fitmaxstab(rain,coord[,1:2],"twhitmat",loc.form,scale.form,shape.form)

TIC(ext_stab,ext_exp,ext_cauchy,ext_matern)

## SLIDE 29 - Diagnostics
plot(ext_stab)

## SLIDE 30 - Additional plots
#f-madogram
par(mar=c(4,5,1,1))
fmadogram(rain, coord, ext_stab, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, ext_stab, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)

#simulating from the estimated copula
sim_text <- rmaxstab(1,as.matrix(coord_sim),cov.mod="powexp",smooth=ext_stab$param[3],nugget=ext_stab$param[1],range=ext_stab$param[2],DoF=ext_stab$param[4])
#define the spatial model for the GEV parameters
param.loc <- ext_stab$param[5] + ext_stab$param[6] * coord_sim[,1]+ ext_stab$param[7]*coord_sim[,2]
param.scale <- ext_stab$param[8] + ext_stab$param[9] * coord_sim[,1]+ ext_stab$param[10]*coord_sim[,2]
param.shape <- ext_stab$param[11]
#creating the data
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(sim_text[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)


# SLIDE 31 - Fitting max-stable models
sch_stab <- fitmaxstab(rain,coord[,1:2],"powexp",loc.form,scale.form,shape.form)
sch_exp <- fitmaxstab(rain,coord[,1:2],"powexp",loc.form,scale.form,shape.form,smooth=1)
sch_cauchy <- fitmaxstab(rain,coord[,1:2],"cauchy",loc.form,scale.form,shape.form)
sch_matern <- fitmaxstab(rain,coord[,1:2],"whitmat",loc.form,scale.form,shape.form)

TIC(sch_stab,sch_exp,sch_cauchy,sch_matern)

brown <- fitmaxstab(rain,coord[,1:2],"brown",loc.form,scale.form,shape.form)

smith <- fitmaxstab(rain,coord[,1:2],"gauss",loc.form,scale.form,shape.form,iso=T)

gg_stab <- fitmaxstab(rain,coord[,1:2],"gpowexp",loc.form,scale.form,shape.form)
gg_exp <- fitmaxstab(rain,coord[,1:2],"gpowexp",loc.form,scale.form,shape.form,smooth=1)
gg_cauchy <- fitmaxstab(rain,coord[,1:2],"gcauchy",loc.form,scale.form,shape.form)
gg_matern <- fitmaxstab(rain,coord[,1:2],"gwhitmat",loc.form,scale.form,shape.form)

TIC(gg_stab,gg_exp,gg_cauchy,gg_matern)

TIC(smith,sch_exp,brown,gg_stab)

## SLIDE 32 - Madograms
par(mar=c(4,5,1,1))
fmadogram(rain, coord, smith, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, smith, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)

par(mar=c(4,5,1,1))
fmadogram(rain, coord, sch_exp, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, sch_exp, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)

par(mar=c(4,5,1,1))
fmadogram(rain, coord, brown, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, brown, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)

par(mar=c(4,5,1,1))
fmadogram(rain, coord, gg_stab, which = "ext",lwd=2,pch=19,cex=0.4,col="gray")
fmadogram(rain, coord, gg_stab, n.bin=300, which = "ext",lwd=2,pch=19,cex=0.4,add=T)


## SLIDE 32 - Simulations
# Smith model
simsmith <- rmaxstab(1,as.matrix(coord_sim),cov.mod="gauss",cov11=smith$param[1],cov12=smith$param[2],cov22=smith$param[3])
param.loc <- smith$param[4] + smith$param[5] * coord_sim[,1]+ smith$param[6]*coord_sim[,2]
param.scale <- smith$param[7] + smith$param[8] * coord_sim[,1]+ smith$param[9]*coord_sim[,2]
param.shape <- smith$param[10]
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(simsmith[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

# Geometric Gaussian
simgg <- rmaxstab(1,as.matrix(coord_sim),cov.mod="gpowexp",sigma2=gg_stab$param[1],nugget=gg_stab$param[2],range=gg_stab$param[3],smooth=gg_stab$param[4])
param.loc <- gg_stab$param[5] + gg_stab$param[6] * coord_sim[,1]+ gg_stab$param[7]*coord_sim[,2]
param.scale <- gg_stab$param[8] + gg_stab$param[9] * coord_sim[,1]+ gg_stab$param[10]*coord_sim[,2]
param.shape <- gg_stab$param[11]
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(simgg[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

#Schlather model
simsch <- rmaxstab(1,as.matrix(coord_sim),cov.mod="powexp",nugget=sch_exp$param[1],range=sch_exp$param[2],smooth=1)
param.loc <- sch_exp$param[4] + sch_exp$param[5] * coord_sim[,1]+ sch_exp$param[6]*coord_sim[,2]
param.scale <- sch_exp$param[7] + sch_exp$param[8] * coord_sim[,1]+ sch_exp$param[9]*coord_sim[,2]
param.shape <- sch_exp$param[10]
data <- matrix(NA, n.obs, n^2)
for (i in 1:n^2)
  data[,i] <- frech2gev(simsch[1,i], param.loc[i], param.scale[i], param.shape)
dat <-matrix(0,n,n)
for(i in 1:n){
  dat[,i]<- data[1,(n*(i-1)+1):(n*i)]
}
image.plot(x,y,dat,xlab="Lon",ylab="Lat")
swiss(add=TRUE,city=TRUE)
points(coord,pch=19,cex=0.6)

## SLIDE 34 - Diagnostics
plot(gg_stab)

## SLIDE 35 - Return levels
SpatialExtremes::map(gg_stab,x,y,param = "loc",xlab="Lon",ylab="Lat")
swiss(add=T,city=T)
points(coord,pch=19,cex=0.6)

SpatialExtremes::map(gg_stab,x,y,param = "scale",xlab="Lon",ylab="Lat")
swiss(add=T,city=T)
points(coord,pch=19,cex=0.6)

SpatialExtremes::map(gg_stab,x,y,ret.per=25,xlab="Lon",ylab="Lat")
swiss(add=T,city=T)
points(coord,pch=19,cex=0.6)

## SLIDE 36 - Conditional return levels
cond.cord <- c(683.26,248)
condmap(gg_stab,cond.cord,x,y)
swiss(add=T,city=T)
points(coord,pch=19,cex=0.6)

## SLIDE 38 - Concurrence probabilities
concprob(rain,coord[,1:2], smith,which="kendall",col="grey",pch=19,cex=0.5)
concprob(rain,coord[,1:2], which="kendall",n.bins=300,pch=19,cex=0.5,add=T)

## SLIDES 39/41 - Conditional Concurrence map
# plot for location 78 - one can try other locations too
cmap <- concurrencemap(rain,coord[,1:2],type=78,n.grid = 500,plot.border = swiss,plot=F) 
image.plot(cmap$x,cmap$y,cmap$z,xlab="lon",ylab="lat")
points(coord[78,1],coord[78,2],pch=4,cex=2)
swiss(add=T,city=T)

## SLIDE 42 - Concurrence map
cmap <- concurrencemap(rain,coord[,1:2],ngrid=1000,plot=F)
image.plot(cmap$x,cmap$y,cmap$z,xlab="lon",ylab="lat")
swiss(add=T,city=T)
points(coord,pch=19,cex=0.6)

