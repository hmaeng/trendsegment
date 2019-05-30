
require(trendsegmentR)
require(IDetect)
require(not)
require(genlasso)
require(xtable)
require(dplyr)
library(tidyr)
require(anomaly)




####################################################################
#   1. Iceland Temperature
#   http://berkeleyearth.org/data/
####################################################################

MyData <- read.csv(file="GlobalLandTemperaturesByCity.csv", header=TRUE, sep=",")
CR <- MyData[MyData$Country=="Iceland",]
CT <- CR[CR$City=="Reykjavík",]
DAT <- CT
d <- data.frame(date = DAT[,1])
DAT <- cbind(separate(d, "date", c("Year", "Month", "Day"), sep = "-"), DAT)
dat <- DAT[DAT$Month=="01" & DAT$Year>=1763, 1:6, drop=F] # January
head(dat)

### Figure 3-(a)
par(mfrow=c(2,2),mar=c(4.5,4,2,2))
plot(dat[,1], dat[,5], type="b",
  col=1, pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
symbols(x=c(1918), y=c(dat[which(dat[,1]==1918),5]), circles=3, col=2, add=T, inches=F)

### Figure 3-(b)
plot(dat[,1], dat[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
tsfit <- ts(x=dat[,5], thr=1.3, p=0.04, bal=0)
c(dat[tsfit$cpts,1])
lines(dat[,1], tsfit$fit, col=1, lty=1, lwd=2) # TS
legend("bottomright", c("obs","TrendSegment"), pch=c(19, rep(NA, 1))
  , col=c("darkgray",1), lty=c(1,1), lwd=c(1,2), cex=1, bty="n", ncol=1)

### Figure 3-(c)
plot(dat[,1], dat[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
lines(dat[,1], not.sic(dat[,5])$fit, col=2, lwd=2, lty=1) # NOT
lines(dat[,1], id(dat[,5])$fit, col=4, lwd=2, lty=2) # ID
legend("bottomright", c("NOT","ID")
  , col=c(2,4), lty=c(1,2), lwd=rep(2,2), cex=1, bty="n", ncol=1)
c(dat[not.sic(dat[,5])$cpts,1])
c(dat[id(dat[,5])$cpts,1])

### Figure 3-(d)
plot(dat[,1], dat[,5], type="b",
  col="darkgray", pch=19, ylab="", xlab="year", cex.axis=1.5, cex.lab=1.5)
lines(dat[,1], tf(dat[,5])$fit, col=3, lwd=2, lty=1) # TF
lines(dat[,1], cpop(dat[,5])$fit, col=2, lty=4, lwd=2) # CPOP
legend("bottomright", c("TF","CPOP")
  , col=c(3:2), lty=c(1,4), lwd=rep(2,2), cex=1, bty="n", ncol=1)
c(dat[tf(dat[,5])$cpts,1])
c(dat[cpop(dat[,5])$cpts,1])





####################################################################
#   2. Sea Ice Extent of Arctic and Antarctic
#   http://nsidc.org/data/nsidc-0051.html
#   https://www.kaggle.com/nsidcorg/daily-sea-ice-extent-data
#   ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/south/daily/data/
#   ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/
####################################################################

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



a <- 9 # month=9 (September) is the coldest in Arctic and hottest in Antarctic
#a <- 2 # month=2 (February) is reverse



### Arctic
north$ym <- paste(north$Year, north$Month)
north.mean <- north %>% group_by(ym) %>% summarise(mean.extent = mean(Extent))
north.mean <- separate(north.mean, ym, into = c('year', 'month'), sep=" ")
north.eachmonth <- north.mean[north.mean$month==a,]
x <- north.eachmonth$mean.extent

### Figures 4-(a),(c) in the paper
plot(north.eachmonth$year, north.eachmonth$mean.extent,
  type="p",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(north.eachmonth$mean.extent)[1]-0.5, range(north.eachmonth$mean.extent)[2]+0.5))

lines(north.eachmonth$year, ts(x, thr=1.3, p=0.04, bal=0)$fit, col=1, lty=1, lwd=2) # TS

### Figures 1-(a) and 2-(a) in the supplementary materials
plot(north.eachmonth$year, north.eachmonth$mean.extent,
  type="p",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(north.eachmonth$mean.extent)[1]-0.5, range(north.eachmonth$mean.extent)[2]+0.5))

lines(north.eachmonth$year, not.sic(x)$fit, col=2, lwd=2, lty=1) # NOT
lines(north.eachmonth$year, id(x)$fit, col=4, lwd=2, lty=2) # ID

legend("topright", c("NOT","ID")
  , col=c(2,4), lty=c(1,2), lwd=rep(2,2), cex=1.5, bty="n", ncol=1)

### Figures 1-(b) and 2-(b) in the supplementary materials
plot(north.eachmonth$year, north.eachmonth$mean.extent,
  type="p",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(north.eachmonth$mean.extent)[1]-0.5, range(north.eachmonth$mean.extent)[2]+0.5))

lines(north.eachmonth$year, tf(x)$fit, col=3, lwd=2, lty=1) # TF
lines(north.eachmonth$year, cpop(x)$fit, col=2, lty=2, lwd=2) # CPOP

legend("topright", c("TF","CPOP")
  , col=c(3:2), lty=c(1,2), lwd=rep(2,2), cex=1.5, bty="n", ncol=1)




### Antarctic
south$ym <- paste(south$Year, south$Month)
south.mean <- south %>% group_by(ym) %>% summarise(mean.extent = mean(Extent))
south.mean <- separate(south.mean, ym, into = c('year', 'month'), sep=" ")
south.eachmonth <- south.mean[south.mean$month==a,]
x <- south.eachmonth$mean.extent

### Figures 4-(b),(d) in the paper
plot(south.eachmonth$year, south.eachmonth$mean.extent,
  type="b",
  xlab="year", ylab="ice extent", col="darkgray", pch=19,cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(south.eachmonth$mean.extent)[1]-0.5, range(south.eachmonth$mean.extent)[2]+0.5))

tsfit <- trendsegment(x=x, p=0.04, th.const = 1.3, bal=0) # TS
lines(south.eachmonth$year, tsfit$est, col=1, lty=1, lwd=2) # TS
south.eachmonth$year[tsfit$cpt]

### Figures 3-(a) and 4-(a) in the supplementary materials
plot(south.eachmonth$year, south.eachmonth$mean.extent,
  type="p",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(south.eachmonth$mean.extent)[1]-0.5, range(south.eachmonth$mean.extent)[2]+0.5))

lines(south.eachmonth$year, not.sic(x)$fit, col=2, lwd=2, lty=1) # NOT
lines(south.eachmonth$year, id(x)$fit, col=4, lwd=2, lty=2) # ID

legend("topleft", c("NOT","ID")
  , col=c(2,4), lty=c(1,2), lwd=rep(2,2), cex=1.5, bty="n", ncol=1)

### Figures 3-(b) and 4-(b) in the supplementary materials
plot(south.eachmonth$year, south.eachmonth$mean.extent,
  type="p",
  xlab="year", ylab="ice extent", col="darkgray", pch=19, cex.axis=1.5, cex.lab=1.5,
  ylim=c(range(south.eachmonth$mean.extent)[1]-0.5, range(south.eachmonth$mean.extent)[2]+0.5))

lines(south.eachmonth$year, tf(x)$fit, col=3, lwd=2, lty=1) # TF
lines(south.eachmonth$year, cpop(x)$fit, col=2, lty=2, lwd=2) # CPOP

legend("topleft", c("TF","CPOP")
  , col=c(3:2), lty=c(1,2), lwd=rep(2,2), cex=1.5, bty="n", ncol=1)



