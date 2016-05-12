#### A WORKFLOW FOR REPRODUCING MEAN BENTHIC GAS FLUXES ####

# Copyright 2016 Robinson W. Fulweiler, Hollie E. Emery, Timothy J. Maguire ####

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
setEPS(width=190/25.4,height=115/25.4)
postscript("Figure2.eps")

#plot A - Providence River Estuary
par(fig=c(.15,.55,.2,.9),mai = c(0,0,0,0))
mean <- N2$ProvN2_mean
sderr <- N2$ProvN2_sderr
names <- N2$datelabel
plotter(mean,sderr,names,"green4",c(-300,300))
axis(2,cex.axis=.7,tcl=-.2,las=2) #adds an axis with labels
mtext(side = 2, text = expression(Net ~ Sediment ~ N[2]-N ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.8, cex=.7)
text(1,250,"a")
text(12,-150, bquote(Mean: .(round(mean(mean,na.rm=T))) %+-% .(round(se(mean))) ~ mu ~ mol ~ m^{-2} ~ h^{-1}), cex=.7, col = "green4", adj=c(0,0))

#plot B - Narragansett Bay
par(fig=c(0.55,.95,.2,.9),mai=c(0,0,0,0),new=TRUE)
mean <- N2$BayN2_mean
sderr <- N2$BayN2_sderr
names <- N2$datelabel
plotter(mean,sderr,names,"red3",c(-300,300))
axis(2,labels=FALSE,tcl=-.2) #adds an axis without labels
text(1,250,"b")
text(12,-150, bquote(Mean: .(round(mean(mean,na.rm=T))) %+-% .(round(se(mean))) ~ mu ~ mol ~ m^{-2} ~ h^{-1}), cex=.7, col = "red3", adj =c(0,0))

dev.off()

#### Figure 3 ####
setEPS(width=190/25.4,height=115/25.4)
postscript("Figure3.eps")

#plot A - Providence River Estuary
par(fig=c(.15,.55,.2,.9),mai = c(0,0,0,0))
mean <- O2$ProvO2_mean
sderr <- O2$ProvO2_sderr
names <- O2$datelabel
plotter(mean,sderr,names,"green4",c(0,8000))
axis(2,cex.axis=.7,tcl=-.2,las=2) #adds an axis with labels
mtext(side = 2, text = expression(Sediment ~ O[2] ~ Demand ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.8, cex=.7)
text(1,7000,"a",adj=c(0,0))
text(12,7000,bquote(Mean: .(round(mean(mean,na.rm=T))) %+-% .(round(se(mean))) ~ mu ~ mol ~ m^{-2} ~ h^{-1}), cex=.7, col = "green4",adj = c(0,0))

#plot B - Narragansett Bay
par(fig=c(0.55,.95,.2,.9),mai=c(0,0,0,0),new=TRUE)
mean <- O2$BayO2_mean
sderr <- O2$BayO2_sderr
names <- O2$datelabel
plotter(mean,sderr,names,"red3",c(0,8000))
axis(2,labels=FALSE,tcl=-.2) #adds an axis without labels
text(1,7000,"b", adj=c(0,0))
text(12,7000,bquote(Mean: .(round(mean(mean,na.rm=T))) %+-% .(round(se(mean))) ~ mu ~ mol ~ m^{-2} ~ h^{-1}), cex=.7, col = "red3",adj = c(0,0))

dev.off()

#### Figure 4 ####

lmo2 = mimsdata$O2Flux/2
prov_o2 = mimsdata[which(mimsdata$Site == "Prov"),]$O2Flux/2
mid_o2 = mimsdata[which(mimsdata$Site == "Mid Bay"),]$O2Flux/2
Prov_fit = lm(mimsdata[which(mimsdata$Site == "Prov"),]$N2Flux~prov_o2)
Mid_fit = lm(mimsdata[which(mimsdata$Site == "Mid Bay"),]$N2Flux~mid_o2)
N2_fit = lm(mimsdata$N2Flux~lmo2)

#make a vector of plotting symbols
symbol=ifelse(mimsdata[,1]=="Prov",19,ifelse(mimsdata[,1]=="Mid Bay",1,NA))

setEPS(width=190/25.4,height=115/25.4)
postscript("Figure4.eps")
par(mai=c(1.02,0.82,0.82,0.42))
plot(lmo2, mimsdata$N2Flux, pch=symbol,  ylim = c(-300,350),xlab="", ylab="")
abline(N2_fit)
abline(Prov_fit,col="green", lty=3)
abline(Mid_fit, col="blue", lty=3)
mtext(side = 2, text = expression(Net ~ Sediment ~ N[2]-N ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 1.8)
mtext(side = 1, text = expression(O[2] ~ Flux ~ (mu ~ mol ~ m^{-2} ~ h^{-1})), line = 2.5)
options(scipen=999)
text(4100,300,bquote(atop(y == .(round(N2_fit$coeff[2],2))*x + .(round(N2_fit$coeff[1],2)), R^{2} ==.(round(summary(N2_fit)$r.squared,2))~~~~ p ==.(round(anova(N2_fit)$'Pr(>F)'[1],4)))), adj=1) 
dev.off()
