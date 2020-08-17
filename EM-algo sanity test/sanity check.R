#------------ Summary and objectives --------------------
#---- //objectives: * explain how EM-algo works
#                   * what happens if we overest./underest.
#---- //@import: PSA_EM algo.R
#---- //@import: package::dplyr
#----------------------------------------------------end
require(dplyr)
rootProject <- "/Users/Vita/sfuvault/RProjects/Pairwise Survival Analysis"
setwd(rootProject)
source("PSA_EM algo with restriction.R")
setwd(paste(rootProject, "[Pairwise Survival] EM algo analysis", sep="/"))

#------------ Sanity check -------------------------------
#---- //objectives: * generate function to draw random sample with mixture model pdf
#                   * this data implies serial interval
#----------------------------------------------------end
rmixGamma <- function(n, alpha, beta, weight) {
  m <- length(weight)
  U <- runif(n)
  rsamp <- numeric(n); omega <- numeric(n)
  
  for(i in 1:m) {
    x <- 0; j <- 1
    while(j<=i) {x <- x+weight[j]; j <- j+1}
    omega[i] <- x
  }
  
  for(i in 1:n) {
    for(j in 1:m) {
      if(U[i]<omega[j]) {rsamp[i] <- rgamma(1, j*alpha, 1/beta)}
    }
  }
  return(rsamp)
}
EMsanityCheck <- function(nrep=50, k=4, discretization=F, rounding=F) {
  nobs <- 100
  result <- NULL
  i <- 1
  
  while(i<=nrep) {
    shape <- runif(1,1,3); scale <- runif(1,2,5)
    weight <- runif(k)
    weight <- weight/sum(weight)
    rsamp <- rmixGamma(nobs, shape, scale, weight)
    if(rounding) rsamp <- round(rsamp)
    if(discretization) {estimates <- EMgamma(rsamp, k); Discrete. <- c(1,1,1)}
    else {estimates <- EMgammaND(rsamp, k); Discrete. <- c(0,0,0)}
    
    if(estimates$niter<1000 & max(estimates$weights)>0.7) {
      record <- data.frame(true_alpha=shape, true_beta=scale,
                           est_alpha=estimates$parameter[1],
                           est_beta=estimates$parameter[2],
                           true_mean=shape*scale,
                           est_mean=prod(estimates$parameter),
                           crude_mean=mean(rsamp), niter=estimates$niter,
                           w1=estimates$weights[1], w2=estimates$weights[2],
                           w3=estimates$weights[3], w4=estimates$weights[4],
                           status=ifelse(discretization, "with discrete.", "without discrete."))
      result <- bind_rows(result,record)
      i <- i+1
    }
    #if((i %% 100)==0) cat("iteration left:", nrep-i, "\n")
  }
  linregAlpha <- lm(est_alpha~true_alpha-1, data=result)
  linregBeta <- lm(est_beta~true_beta-1, data=result)
  linregMean <- lm(est_mean~true_mean-1, data=result)
  param <- rbind(summary(linregAlpha)$coef, summary(linregBeta)$coef, summary(linregMean)$coef)
  param <- cbind(param, Discrete.)
  
  out <- list(record=result, summary=param)
  class(out) <- "EM-check"
  return(out)
}

#------------ Over/Under-estimation check --------------------------
#---- //objectives: * what happens if we estimate with incorrect k
#---------------------------------------------------------------end
EMoverUnderEst <- function(nrep=50, k=1, discretization=F, rounding=F) {
  nobs <- 100
  result <- NULL
  
  for(i in 1:nrep) {
    shape <- runif(1,1,5); scale <- runif(1,1,5)
    weight <- runif(k); weight <- weight/sum(weight)
    rsamp <- rmixGamma(nobs, shape, scale, weight)
    if(rounding) rsamp <- round(rsamp)
    
    if(discretization) {
      for(j in 1:4) {
        est <- EMgamma(rsamp, j)
        record <- data.frame(true_alpha=shape, true_beta=scale,
                             est_alpha=est$parameter[1], est_beta=est$parameter[2],
                             true_mean=shape*scale, 
                             est_mean=prod(est$parameter),
                             crude_mean=mean(rsamp), 
                             w1=est$weights[1], w2=est$weights[2], w3=est$weights[3], w4=est$weights[4],
                             status=ifelse(discretization, "with discrete.", "without discrete."),
                             k=paste("ncomp: k=",k, sep=""), EMcomp=paste("EMcomp=",j, sep=""))
        result <- bind_rows(result, record)
      }
    } else {
      for(j in 1:4) {
        est <- EMgammaND(rsamp, j)
        record <- data.frame(true_alpha=shape, true_beta=scale,
                             est_alpha=est$parameter[1], est_beta=est$parameter[2],
                             true_mean=shape*scale, 
                             est_mean=prod(est$parameter),
                             crude_mean=mean(rsamp), 
                             w1=est$weights[1], w2=est$weights[2], w3=est$weights[3], w4=est$weights[4],
                             status=ifelse(discretization, "with discrete.", "without discrete."),
                             k=paste("ncomp: k=",k, sep=""), EMcomp=paste("EMcomp=",j, sep=""))
        result <- bind_rows(result, record)
      }
    }
  }
  linregAlpha1 <- lm(true_alpha~est_alpha-1, data=(result %>% filter(EMcomp=="EMcomp=1")))
  linregAlpha2 <- lm(true_alpha~est_alpha-1, data=(result %>% filter(EMcomp=="EMcomp=2")))
  linregAlpha3 <- lm(true_alpha~est_alpha-1, data=(result %>% filter(EMcomp=="EMcomp=3")))
  linregAlpha4 <- lm(true_alpha~est_alpha-1, data=(result %>% filter(EMcomp=="EMcomp=4")))

  linregBeta1 <- lm(true_beta~est_beta-1, data=(result %>% filter(EMcomp=="EMcomp=1")))
  linregBeta2 <- lm(true_beta~est_beta-1, data=(result %>% filter(EMcomp=="EMcomp=2")))
  linregBeta3 <- lm(true_beta~est_beta-1, data=(result %>% filter(EMcomp=="EMcomp=3")))
  linregBeta4 <- lm(true_beta~est_beta-1, data=(result %>% filter(EMcomp=="EMcomp=4")))

  linregMean1 <- lm(true_mean~est_mean-1, data=(result %>% filter(EMcomp=="EMcomp=1", !(is.nan(est_mean)))))
  linregMean2 <- lm(true_mean~est_mean-1, data=(result %>% filter(EMcomp=="EMcomp=2", !(is.nan(est_mean)))))
  linregMean3 <- lm(true_mean~est_mean-1, data=(result %>% filter(EMcomp=="EMcomp=3", !(is.nan(est_mean)))))
  linregMean4 <- lm(true_mean~est_mean-1, data=(result %>% filter(EMcomp=="EMcomp=4", !(is.nan(est_mean)))))

  temp <- cbind(k=rep(k,12), ncomp=rep(c(1:4), times=3))
  param <- rbind(summary(linregAlpha1)$coef, summary(linregAlpha2)$coef,
                 summary(linregAlpha3)$coef, summary(linregAlpha4)$coef,
                 summary(linregBeta1)$coef, summary(linregBeta2)$coef,
                 summary(linregBeta3)$coef, summary(linregBeta4)$coef,
                 summary(linregMean1)$coef, summary(linregMean2)$coef,
                 summary(linregMean3)$coef, summary(linregMean4)$coef)
  param <- cbind(temp, param)
  
  out <- list(record=result, summary=param)
  class(out) <- "EM-check"
  return(out)
}


