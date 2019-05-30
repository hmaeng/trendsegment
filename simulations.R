
require(trendsegmentR)
require(IDetect)
require(not)
require(genlasso)
require(xtable)
require(dplyr)
# CPOP R code can be imported from:
# http://www.research.lancs.ac.uk/portal/en/datasets/cpop(56c07868-3fe9-4016-ad99-54439ec03b6c).html
source("CPOP.R")
source("estimateCPOP.R")


have.signal <- function(model){
  if(model$cpt.type == "pcwsConstMean"){

    signal <- rep(0, model$n)

    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    signal[segments[1,1]:segments[1,2]] <- model$start[1]

    for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]

  }else if(model$cpt.type == "pcwsLinContMean" & length(model$cpt)>0){
    signal <- rep(0, model$n)
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))

    slope <- model$start[2]
    signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
    for(j in 2:nrow(segments)) {

      slope <- slope +  model$jump.size[j-1]

      for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
    }

  }else if(model$cpt.type == "pcwsLinContMean" & length(model$cpt)==0){
    signal <- rep(0, model$n)
    segments <- matrix(c(1, model$n), ncol=2)
    slope <- model$start[2]
    signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]

  }else if(model$cpt.type == "pcwsLinMean"){

    signal <- rep(0, model$n)
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))

    slope <- model$start[2]
    signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope

    for(j in 2:nrow(segments)) {

      slope <- slope +  model$jump.size[j-1,2]

      signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
      if(segments[j,1]+1 < segments[j,2]){
        for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
      }
    }
  }

  return(signal)
}

finding.dH <- function(chp, modelnum, models){
  n <- models[[modelnum]]$n
  segments.endpoints.est <- sort(unique(c(0, chp, n)))
  segments.endpoints.true <- sort(unique(c(0, models[[modelnum]]$cpt, n)))

  distm <- abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)), nrow=length(segments.endpoints.est))
    -matrix(rep(segments.endpoints.true, length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=T))

  screening.dist <- max(apply(distm, 2, min)) * 100 / n
  precision.dist <- max(apply(distm, 1, min)) * 100 / n
  dH <- mean((abs(screening.dist-precision.dist) + screening.dist + precision.dist)/2)

  return(dH)
}

assess <- function(object, modelnum, models){
  qdiff <- length(which(object$cpts>0)) - length(models[[modelnum]]$cpt)
  mse <- mean((have.signal(models[[modelnum]])-object$fit)^2)
  dh <- finding.dH(chp=object$cpts, modelnum=modelnum, models=models)
  return(list(qdiff=qdiff, mse=mse, dh=dh))
}


#######################################
#   Models
#######################################

model.wave1 <-  list(name = "wave1",
  cpt.type = "pcwsLinContMean",
  cpt = (1:9) * 150,
  jump.size = (-1)^{1:16} / 25,
  n = 150* 10,
  start=c(-1,1/50))

model.wave2 <-  list(name = "wave2",
  cpt.type = "pcwsLinMean",
  cpt = (1:20) * 60,
  jump.size = matrix(c(rep(c(1, -1),10), (-1)^{1:20}/8), ncol=2),
  n = 20*60+60,
  start=c(-1,1/16))

model.mix1 <-  list(name = "mix1",
  cpt.type = "pcwsLinMean",
  cpt = (1:7) * 256,
  jump.size = matrix(c(c(0,-1,0,-1,1,1,0), c(1,-1,-1,1,1,-2,2)/64), ncol=2),
  n = 2048,
  start=c(0,0))

model.mix2 <-  list(name = "mix2",
  cpt.type = "pcwsLinMean",
  cpt = sort(c((1:7)*256, 257, 1793)),
  jump.size = matrix(c(c(-7,7,2,-2,1,-1,1,-7,7), c(0,-1,1,1,-2,2.5,-1.5,0,-1.5)/64), ncol=2),
  n = 2048,
  start=c(2,0))

model.mix3 <- list(name = "mix3",
  cpt.type = "pcwsLinMean",
  cpt = sort(c((1:7)*256, 512+30, 1280+30, 1793)),
  jump.size = matrix(c(c(-4, 6,-4,-1,1,-5,4,2,-7,7), c(0, 2,-3,2,-2,1,0,1,0,-2)/64), ncol=2),
  n = 2048,
  start=c(2,0))

model.linsgmts <-  list(name = "linsgmts", # bump function length=10 / jump size=2
  cpt.type = "pcwsLinMean",
  cpt = sort(c((1:4)*2*256, (1:4)*2*256+5)),
  jump.size = matrix(c(rep(c(6,-6-4/64), 4), rep(c(1,-1)/64,4)), ncol=2),
  n = 256*9,
  start=c(-1,0))

model.teeth <-  list(name = "teeth",
  cpt.type = "pcwsConstMean",
  cpt = (1:7) * 100,
  jump.size = 2*(-1)^{1:7},
  n = 100 * 8,
  start = 1)

model.lin <-  list(name = "lin",
  cpt.type = "pcwsLinContMean",
  cpt = c(),
  jump.size = c(),
  n = 150* 10,
  start=c(-1, 2/1500))

models <- list(model.wave1, model.wave2 , model.mix1,
  model.mix2, model.mix3,
  model.linsgmts, model.teeth, model.lin)

par(mfrow=c(4,2),mar=rep(3,4))

for(i in 1:length(models)){
  truex <- have.signal(models[[i]])
  x <- truex + rnorm(length(truex), sd=1)
  plot(x, type="l", ylab="", xlab="t", ylim=range(c(x, truex)), col=8, main=paste("( M",i,")" ,models[[i]]$name))
  lines(truex, col=1)
}

#######################################
#   Simulations
#######################################
sims<- function(N, modelnum, sgma, thr, p, bal, noise){

  truex <- have.signal(models[[modelnum]]) # signal only
  n <- length(truex)

  par(mfrow=c(1,1),mar=c(2,2,1.5,1.5))

  ### Assesement
  result <- list(qdiff=matrix(NA, nrow=N, ncol=6), mse=matrix(NA, nrow=N, ncol=6), dh=matrix(NA, nrow=N, ncol=6), time=matrix(NA, nrow=N, ncol=6))
  result <- lapply(result, function(x) {colnames(x) <- c("TS", "NOT", "ID", "TF", "CPOP", "BUP"); x})

  for(K in 1:N){
    set.seed(K)

    if(noise=="norm1"){
      x <- truex + rnorm(n, sd=sgma)
    } else if(noise=="t5"){
      x <- truex + rt(n, df=5) * sqrt(3/5)
    } else if(noise=="laplace"){
      x <- truex + extraDistr::rlaplace(n, mu = 0, sigma = 1)
    } else if(noise=="ar103"){
      x <- truex + arima.sim(n=n, list(ar=c(0.3)), rand.gen = function(n) rnorm(n, sd=sgma))
    }

    plot(x, type="l", col=8, xlab="t", ylab="",
      ylim=range(x),
      main=paste(models[[modelnum]]$name, " ( seed = ", K, " )"))
    lines(truex, col=7, lwd=2)


    #####################################################################
    ################################## TS ###############################
    #####################################################################
    obj <- ts(x, thr=thr, p=p, bal=bal) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,1] <- a.obj$qdiff
    result[[2]][K,1] <- a.obj$mse
    result[[3]][K,1] <- a.obj$dh
    result[[4]][K,1] <- obj$elapsed

    lines(obj$fit, col=1, lwd=1, lty=1)
    #####################################################################
    ################################# NOT ###############################
    #####################################################################
    obj <- not.sic(x) ## cpts is turned "NA" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,2] <- a.obj$qdiff
    result[[2]][K,2] <- a.obj$mse
    result[[3]][K,2] <- a.obj$dh
    result[[4]][K,2] <- obj$elapsed

    lines(obj$fit, col=2, lwd=1, lty=2)
    #####################################################################
    ################################# ID ################################
    #####################################################################
    obj <- id(x) ## cpts is turned "0" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,3] <- a.obj$qdiff
    result[[2]][K,3] <- a.obj$mse
    result[[3]][K,3] <- a.obj$dh
    result[[4]][K,3] <- obj$elapsed

    lines(obj$fit, col=3, lwd=1, lty=3)
    #####################################################################
    ########################## Trend Filtering ##########################
    #####################################################################
    obj <- tf(x) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,4] <- a.obj$qdiff
    result[[2]][K,4] <- a.obj$mse
    result[[3]][K,4] <- a.obj$dh
    result[[4]][K,4] <- obj$elapsed

    lines(obj$fit, col=4, lwd=1, lty=4)
    #####################################################################
    ################################ CPOP ###############################
    #####################################################################
    obj <- cpop(x) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,5] <- a.obj$qdiff
    result[[2]][K,5] <- a.obj$mse
    result[[3]][K,5] <- a.obj$dh
    result[[4]][K,5] <- obj$elapsed

    lines(obj$fit, col=5, lwd=1, lty=5)

    #####################################################################
    ################################# BUP ###############################
    #####################################################################
    obj <- bup(x) ## cpts is turned "numeric(0)" if nothing is detected
    a.obj <- assess(object=obj, modelnum=modelnum, models=models)
    result[[1]][K,6] <- a.obj$qdiff
    result[[2]][K,6] <- a.obj$mse
    result[[3]][K,6] <- a.obj$dh
    result[[4]][K,6] <- obj$elapsed

    lines(obj$fit, col=6, lwd=1, lty=6)

    legend("topleft", c("obs","true","TS","NOT","ID","TF","CPOP", "BUP"), col=c(8,7,1:6), lty=c(1,1,1:6), lwd=c(1,2,rep(1,6)), cex=0.7, bty="n", ncol=2)
  }

  return(result)

}

m1n1 <- sims(N=100, modelnum=1, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m2n1 <- sims(N=100, modelnum=2, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m3n1 <- sims(N=100, modelnum=3, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m4n1 <- sims(N=100, modelnum=4, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m5n1 <- sims(N=100, modelnum=5, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m6n1 <- sims(N=100, modelnum=6, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m7n1 <- sims(N=100, modelnum=7, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")
m8n1 <- sims(N=100, modelnum=8, sgma=1, thr=1.3, p=0.04, bal=0, noise="norm1")

