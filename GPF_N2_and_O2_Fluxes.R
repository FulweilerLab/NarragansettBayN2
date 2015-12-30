#### Set up ####
## packages that are needed for this script
library(plotrix)
library(plyr)
library(repmis)

## functions to define

# to find standard error
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# to make barplots
plotter <- function(mean, sderr, names, col, ylim){
  plotCI(barplot(mean,col=col,ylim=ylim,yaxt="n",names=names,las=2,cex.names=.7), #makes a barplot of the means with no y axis
  mean,sderr,pch=NA,sfrac=0,add=TRUE) #adds the error bars
  abline(h=0)
  abline(h=mean(mean,na.rm=TRUE),lty=2,lwd=1.5)
}

## Downloads and reads the data off a public figshare link
mimsdata <- repmis::source_data("https://ndownloader.figshare.com/files/3639429?private_link=dbb198fc42b811e58ba406ec4b8d1f61",
                                sep = ",",
                                header = TRUE)

#### Data processing ####

### fix the dates, then redisplay them as year-month format for easy reading
mimsdata$date <- as.POSIXct(strptime(mimsdata$Date,format="%m/%d/%Y",tz="EST"))
mimsdata$yearmonth <- as.Date(paste(format(mimsdata$date,"%Y-%m"),"-01",sep="")) #this makes them all appear on the first of the month to easy aggregating

### find the flux means and standard errors per sample date at each site
## N2
N2_mean <- adply(tapply(mimsdata$N2Flux,list(mimsdata$yearmonth,mimsdata$Site),mean),1)
N2_sderr <- adply(tapply(mimsdata$N2Flux,list(mimsdata$yearmonth,mimsdata$Site),se),1)

## O2
O2_mean <- adply(tapply(mimsdata$O2Flux,list(mimsdata$yearmonth,mimsdata$Site),mean),1)
O2_sderr <- adply(tapply(mimsdata$O2Flux,list(mimsdata$yearmonth,mimsdata$Site),se),1)

# final dataframes for plotting
N2 <- merge(N2_mean,N2_sderr,by="X1",all=TRUE)
colnames(N2) <- c("date","BayN2_mean","ProvN2_mean","BayN2_sderr","ProvN2_sderr")
N2$datelabel <- format(as.Date(N2$date),format="%b-%y")

O2 <- merge(O2_mean,O2_sderr,by="X1",all=TRUE)
colnames(O2) <- c("date","BayO2_mean","ProvO2_mean","BayO2_sderr","ProvO2_sderr")
O2$datelabel <- format(as.Date(O2$date),format="%b-%y")

#### Figure 2 ####

par(mfrow = c(1,2))

#plot A - Providence River Estuary
par(mai = c(1.4,1.4,1,0.1))
mean <- N2$ProvN2_mean
sderr <- N2$ProvN2_sderr
names <- N2$datelabel
plotter(mean,sderr,names,"green4",c(-300,300))
axis(2,cex.axis=.7,tcl=-.2,las=2) #adds an axis with labels
mtext(side = 2, text = expression(Net ~ Sediment ~ N[2]-N ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.8, cex=.7)
text(1,250,"a")

#plot B - Narragansett Bay
par(mai=c(1.4,0.1,1,1.4))
mean <- N2$BayN2_mean
sderr <- N2$BayN2_sderr
names <- N2$datelabel
plotter(mean,sderr,names,"red3",c(-300,300))
axis(2,labels=FALSE,tcl=-.2) #adds an axis without labels
text(1,250,"b")

par(mfrow = c(1,1))


#### Figure 3 ####

#This is to make a similar figure to the above but with O2 
par(mfrow = c(1,2))

#plot A - Providence River Estuary
par(mai=c(1.4,1.4,1,0.1))
mean <- O2$ProvO2_mean
sderr <- O2$ProvO2_sderr
names <- O2$datelabel
plotter(mean,sderr,names,"green4",c(0,8000))
axis(2,cex.axis=.7,tcl=-.2,las=2) #adds an axis with labels
mtext(side = 2, text = expression(Sediment ~ O[2] ~ Demand ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.8, cex=.7)
text(1,7000,"a")

#plot B - Narragansett Bay
par(mai=c(1.4,0.1,1,1.4))
mean <- O2$BayO2_mean
sderr <- O2$BayO2_sderr
names <- O2$datelabel
plotter(mean,sderr,names,"red3",c(0,8000))
axis(2,labels=FALSE,tcl=-.2) #adds an axis without labels
text(1,7000,"b")

par(mfrow = c(1,1))


#### Figure 4 ####

lmo2 = mimsdata$O2Flux/2
N2_fit = lm(mimsdata$N2Flux~lmo2)
#summary(N2_fit)

#make a vector of plotting symbols
symbol=ifelse(mimsdata[,1]=="Prov",19,ifelse(mimsdata[,1]=="Mid Bay",1,NA))

par(mai=c(1.02,0.82,0.82,0.42))
plot(lmo2, mimsdata$N2Flux, pch=symbol,  ylim = c(-300,350),xlab="", ylab="")
abline(N2_fit)
mtext(side = 2, text = expression(Net ~ Sediment ~ N[2]-N ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 1.8)
mtext(side = 1, text = expression(O[2] ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.5)
