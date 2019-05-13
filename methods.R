
#######################################
#   Methods
#######################################

ts <- function(x, thr, bal, method, beta, ...){

  tic <- proc.time()
  object <- trendsegment(x, th.const = thr, bal = bal)
  toc <- proc.time()

  list(fit = object$est, cpts=object$cpt, elapsed=(toc-tic)[3])

}

not.sic <- function(x){

  tic <- proc.time()
  object <- not(x, method = "not", contrast = "pcwsLinMean", parallel = FALSE) # give an error in noiseless input
  cpts <- features(object, penalty="sic")$cpt
  fit <- predict(object, cpt=cpts)
  toc <- proc.time()

  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}

id <- function(x){ # adjusted cpts after fitting
  # e.g. x <- c(0:9, 5, 10:0)+rnorm(22)/3  gives three consecutive change-points
  tic <- proc.time()
  #object <- ID(x, contrast ="slope", ht=F) # give a weird cpts in noiseless input
  object <- ID(x, contrast ="slope", ht=F, lambda=1) # this does not make any difference compared to the upper one in MODEL 5-7 which contain spikes
  if(length(object$cpt)==1 && object$cpt==c(0)){
    cpts <- c()
  } else{
    cpts <- object$cpt
  }

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()

  list(fit = object$fit, cpts=cpts, elapsed=(toc-tic)[3])

}

tf <- function(x){ # adjusted cpts after fitting
  # e.g. x <- c(rep(0,9), 2, 10:0) gives three consecutive cpts
  tic <- proc.time()
  object <-  trendfilter(y=x, ord=1)

  tf.cv <- cv.trendfilter(object)
  cpts <- which(abs(diff(object$fit[,tf.cv$i.min], differences=2)) > sqrt(.Machine$double.eps))+1

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()

  list(fit = object$fit[,tf.cv$i.min], cpts=cpts, elapsed=(toc-tic)[3])

}

cpop <- function(x){ # adjusted cpts after fitting
  # e.g. x <- c(0:9, 5, 10:0)+rnorm(22)/3 returns three consecutive cpts
  tic <- proc.time()
  sig.hat <- mad(diff(diff(x)))/sqrt(6) # give an error in noiseless input
  object <- CPOP.run(x/sig.hat)
  cpts <- object$cpt
  fit <- object$f*sig.hat

  ### adjust change-points as it gives one more change-point for each estimated change-points by finding the point which gives nonzero second derivative.
  if(length(cpts)>2){
    # 1) searching spikes only
    threeeach <- cbind(cpts[1:(length(cpts)-2)], cpts[2:(length(cpts)-1)], cpts[3:(length(cpts))])
    dthr <- t(apply(threeeach,1,diff))
    spike <- which(dthr[,1]==1 & dthr[,2]==1)
    # 2) searching jump only
    twoeach <- cbind(cpts[1:(length(cpts)-1)], cpts[2:(length(cpts))])
    jump <- which(t(apply(twoeach,1,diff))==1)
    jump <-jump[!jump %in% c(spike, spike+1)] # get rid of the cpts for spikes
    # removing redundunt
    if(length(spike)>0 | length(jump)>0){
      cpts <- cpts[-unique(c(spike+2, jump+1))]
    } else{
      cpts <- cpts
    }
  } else if(length(cpts)==2 && diff(cpts)[1]==1){
    cpts <- cpts[-2]
  } else {
    cpts <- cpts
  }

  toc <- proc.time()

  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}

bup <- function(x, max.err=300){

  tic <- proc.time()

  left_x <- seq(1, length(x)-1, by=2)
  right_x <- left_x + 1
  right_x[length(right_x)] <- length(x)
  number_of_segments <- length(left_x)

  segment <- cbind(left_x, right_x, Inf)

  for(i in 1:(number_of_segments-1)){
    lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
    segment[i, 3] <- sum((lmfit$residuals)^2) # sse
  }

  while(min(segment[,3]) < max.err){ # max.err

    i <- which.min(segment[,3])

    if(i==1){
      lmfit <- lm(x[c(segment[i, 1]:segment[i+2, 2])]~ c(segment[i, 1]:segment[i+2, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse

      segment[i, 2] <- segment[i+1, 2]
      segment <- segment[-(i+1),,drop=F]
    } else if(i>1 && i < dim(segment)[1]-1){
      lmfit <- lm(x[c(segment[i, 1]:segment[i+2, 2])]~ c(segment[i, 1]:segment[i+2, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse

      segment[i, 2] <- segment[i+1, 2]
      segment <- segment[-(i+1),,drop=F]

      i <- i-1
      lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse
    } else{
      segment[i, 2] <- segment[i+1, 2]
      segment[i, 3] <- Inf
      segment <- segment[-(i+1),,drop=F]

      i <- i-1
      lmfit <- lm(x[c(segment[i, 1]:segment[i+1, 2])]~ c(segment[i, 1]:segment[i+1, 2]))
      segment[i, 3] <- sum((lmfit$residuals)^2) # sse
    }

  }

  ### change-points
  cpt <- c(segment[-dim(segment)[1], 2])
  ### estimated curve
  est <- rep(NA, length(x))
  for(i in 1:dim(segment)[1]){
    lmfit <- lm(x[segment[i,1]:segment[i,2]]~c(segment[i,1]:segment[i,2]))
    est[c(segment[i,1]:segment[i,2])] <- lmfit$fitted.values
  }


  toc <- proc.time()

  list(fit=est, cpts=cpt, elapsed=(toc-tic)[3])

}
