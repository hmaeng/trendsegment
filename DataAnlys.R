
require(trendsegmentR)
require(IDetect)
require(not)
require(genlasso)
require(xtable)
require(dplyr)
library(tidyr)
require(anomaly)




###############################################
#   1. Iceland Temperature
#   http://berkeleyearth.org/data/
###############################################

MyData <- read.csv(file="GlobalLandTemperaturesByCity.csv", header=TRUE, sep=",")
CR <- MyData[MyData$Country=="Iceland",]
CT <- CR[CR$City=="Reykjavík",]
DAT <- CT
d <- data.frame(date = DAT[,1])
DAT <- cbind(separate(d, "date", c("Year", "Month", "Day"), sep = "-"), DAT)
month <- DAT[DAT$Month=="01", 1:6, drop=F] # January
month <- month[-c(1:max(which(is.na(month[,5])))),,drop=F]
head(month)

##### Figure 3
par(mfrow=c(2,2),mar=c(4.5,4,2,2))
plot(month[,1], month[,5], type="b",
  col=1, pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
symbols(x=c(1918), y=c(month[which(month[,1]==1918),5]), circles=3, col=2, add=T, inches=F)

plot(month[,1], month[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
tsfit <- ts(x=month[,5], thr=1, bal=0)
c(month[tsfit$cpts,1])
lines(month[,1], tsfit$fit, col=1, lty=1, lwd=2) # ts
legend("bottomright", c("obs","TrendSegment"), pch=c(19, rep(NA, 1))
  , col=c("darkgray",1), lty=c(1,1), lwd=c(1,2), cex=1, bty="n", ncol=1)

plot(month[,1], month[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
lines(month[,1], not.sic(month[,5])$fit, col=2, lwd=2, lty=1) # NOT
lines(month[,1], id(month[,5])$fit, col=4, lwd=2, lty=2) # ID
legend("bottomright", c("NOT","ID")
  , col=c(2,4), lty=c(1,2), lwd=rep(2,2), cex=1, bty="n", ncol=1)
c(month[not.sic(month[,5])$cpts,1])
c(month[id(month[,5])$cpts,1])

plot(month[,1], month[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
lines(month[,1], tf(month[,5])$fit, col=3, lwd=2, lty=1) # TF
lines(month[,1], cpop(month[,5])$fit, col=2, lty=4, lwd=2) # CPOP
legend("bottomright", c("TF","CPOP")
  , col=c(3:2), lty=c(1,4), lwd=rep(2,2), cex=1, bty="n", ncol=1)
c(month[tf(month[,5])$cpts,1])
c(month[cpop(month[,5])$cpts,1])





###############################################
#   2. Sea Ice Extent of Arctic and Antarctic
#   http://nsidc.org/data/nsidc-0051.html
#   https://www.kaggle.com/nsidcorg/daily-sea-ice-extent-data
#   ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/south/daily/data/
#   ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/
###############################################

### north
dat <- read.csv(file="N_seaice_extent_daily.csv", header=TRUE, sep=",")
dat <- dat[-1,-6, drop=F] # remove the website address column
# change the factor to numeric
dat$Year <- as.numeric(as.character(dat$Year))
dat$Month <- as.numeric(as.character(dat$Month))
dat$Day <- as.numeric(as.character(dat$Day))
dat$Extent <- as.numeric(as.character(dat$Extent))
north <- dat

### south
dat <- read.csv(file="S_seaice_extent_daily.csv", header=TRUE, sep=",")
dat <- dat[-1,-6, drop=F] # remove the website address column
# change the factor to numeric
dat$Year <- as.numeric(as.character(dat$Year))
dat$Month <- as.numeric(as.character(dat$Month))
dat$Day <- as.numeric(as.character(dat$Day))
dat$Extent <- as.numeric(as.character(dat$Extent))
south <- dat


##### Figure 4
a <- 9 # month=9 (September) is the hottest in Arctic and coldest in Antarctic
#a <- 2 # month=2 (February) is reverse


### Arctic
north$ym <- paste(north$Year, north$Month)
north.mean <- north %>% group_by(ym) %>% summarise(mean.extent = mean(Extent))
north.mean <- separate(north.mean, ym, into = c('year', 'month'), sep=" ")
north.eachmonth <- north.mean[north.mean$month==a,]

par(mfrow=c(1,1),mar=c(5,5,2,2))
plot(north.eachmonth$year, north.eachmonth$mean.extent,
  main=paste("Arctic ( Month = ", a, ")"),
  type="b",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(north.eachmonth$mean.extent)[1]-0.5, range(north.eachmonth$mean.extent)[2]+0.5))
x <- north.eachmonth$mean.extent

lines(north.eachmonth$year, ts(x, thr=1.3, bal=0)$fit, col=1, lty=1, lwd=3) # TS
lines(north.eachmonth$year, not.sic(x)$fit, col=2, lwd=2, lty=2) # NOT
lines(north.eachmonth$year, id(x)$fit, col=3, lwd=2, lty=3) # ID
lines(north.eachmonth$year, tf(x)$fit, col=4, lwd=1, lty=4) # TF
lines(north.eachmonth$year, cpop(x)$fit, col=5, lty=5, lwd=1) # CPOP
legend("topright", c("obs","TS","NOT","ID","TF","CPOP"), pch=c(19, rep(NA, 5))
  , col=c("darkgray",1:5), lty=c(1,1:5), lwd=c(1,rep(2,5)), cex=1, bty="n", ncol=2)


### Antarctic
south$ym <- paste(south$Year, south$Month)
south.mean <- south %>% group_by(ym) %>% summarise(mean.extent = mean(Extent))
south.mean <- separate(south.mean, ym, into = c('year', 'month'), sep=" ")
south.eachmonth <- south.mean[south.mean$month==a,]

par(mfrow=c(1,1),mar=c(5,5,2,2))
plot(south.eachmonth$year, south.eachmonth$mean.extent,
  main=paste("Antarctic ( Month = ", a, ")"),
  type="b",
  xlab="year", ylab="ice extent", col="darkgray", pch=19,cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(south.eachmonth$mean.extent)[1]-0.5, range(south.eachmonth$mean.extent)[2]+0.5))
x <- south.eachmonth$mean.extent

lines(south.eachmonth$year, ts(x, thr=1.3, bal=0)$fit, col=1, lty=1, lwd=3) # ts
lines(south.eachmonth$year, not.sic(x)$fit, col=2, lwd=2, lty=2) # NOT
lines(south.eachmonth$year, id(x)$fit, col=3, lwd=2, lty=3) # ID
lines(south.eachmonth$year, tf(x)$fit, col=4, lwd=1, lty=4) # TF
lines(south.eachmonth$year, cpop(x)$fit, col=5, lty=5, lwd=1) # CPOP
legend("topleft", c("obs","TS","NOT","ID","TF","CPOP"), pch=c(19, rep(NA, 6))
  ,col=c("darkgray",1:5), lty=c(1,1:5), lwd=c(1,rep(2,5)), cex=1, bty="n", ncol=2)
