
library(ks)

library(mgcv)

#####################################################
#####################################################
#####################################################
###########" Simulation 1 Standard Normal ##########
#####################################################
#####################################################
#####################################################



w=150

#####################################################
##############' Model ONE: ReSampling ###############
#####################################################

output1.Sim1 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps = rnorm(n)
  beta0=1
  beta1=3
  y = beta0+beta1*x+eps
  
  L1 <- lm(y ~ x) 
  summary(L1)
  b0 <- L1$coeff[1]
  b1 <- L1$coeff[2]
  
  m <- 2000
  out <- replicate(m, expr={
    i <- sample(1:n, replace=TRUE, size=n)
    xstar <- x[i]
    ystar <- y[i]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
  })
  
  bootCases <- t(out)
  meanCases <- colMeans(bootCases)
  sdCases <- apply(bootCases, 2, "sd")
  
  slope.boot.replicates<- bootCases[,2]
  
  biasInt <- mean(bootCases[,1] - b0)
  biasSlope <- mean(bootCases[,2] - b1)
  rbind(estimate=c(b0, b1), bias=c(biasInt, biasSlope), se=sdCases[1:2], cv=c(biasInt, cv=biasSlope)/sdCases[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.resample.percentile <- quantile(slope.boot.replicates,probs = p.lower,names = F)
  
  ub.resample.percentile <- quantile(slope.boot.replicates,probs = p.upper,names = F)
  
  LB <- as.vector(lb.resample.percentile)
  
  output1.Sim1 <- data.frame(LB)
  
  UB <- as.vector(ub.resample.percentile)
  
  output1.Sim1$UB <- UB
  
  output1.Sim1
  
})


x1.sim1 <- unlist(output1.Sim1)

LB <- x1.sim1[seq(1, length(x1.sim1), 2)]
UB <- x1.sim1[seq(2, length(x1.sim1), 2)]

(C.I_For_ReSampling_Cases_BootStrap.Sim1 <- cbind(LB,UB))

(In_C.I_ReSamp.Sim1 <- sum(C.I_For_ReSampling_Cases_BootStrap.Sim1[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim1[,"UB"]>2.99999999))

(ReSampling.prop.sim1 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim1[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim1[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSampling.Cases.Sim1 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim1[,"UB"]-C.I_For_ReSampling_Cases_BootStrap.Sim1[,"LB"]))


#####################################################
#########' Model Two: Parametric Bootstrap ##########
#####################################################


output2.Sim1 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps = rnorm(n)
  beta0=1
  beta1=3
  y.2 = beta0+beta1*x+eps
  L2 <- lm(y.2 ~ x)
  summary(L2)
  b0.2 <- L2$coeff[1]
  b1.2 <- L2$coeff[2]
  s <- summary(L2)$sigma
  yhat <- L2$fitted.values
  m <- 2000
  
  out.2 <- replicate(m, expr={
    eps.star.2 <- rnorm(n, 0, sd=s)
    xstar.2 <- x
    ystar.2 <- yhat+ eps.star.2
    Lb.2 <- lm(ystar.2 ~ xstar.2)
    s.2 <- summary(Lb.2)$sigma
    c(Lb.2$coeff[1], slope=Lb.2$coeff[2], s=s.2)
  })
  
  bootCases2 <- t(out.2)
  meanCases2 <- colMeans(bootCases2)
  sdCases2 <- apply(bootCases2, 2, "sd")
  
  slope.boot.replicates2<- bootCases2[,2]
  
  biasInt2 <- mean(bootCases2[,1] - b0.2)
  biasSlope2 <- mean(bootCases2[,2] - b1.2)
  rbind(estimate=c(b0.2, b1.2), bias=c(biasInt2, biasSlope2), se=sdCases2[1:2], cv=c(biasInt2, cv=biasSlope2)/sdCases2[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.parametric.percentile <- quantile(slope.boot.replicates2,probs = p.lower,names = F)
  
  ub.parametric.percentile <- quantile(slope.boot.replicates2,probs = p.upper,names = F)
  
  LB <- as.vector(lb.parametric.percentile)
  
  output2.Sim1 <- data.frame(LB)
  
  UB <- as.vector(ub.parametric.percentile)
  
  output2.Sim1$UB <- UB
  
  output2.Sim1
  
})


x2.sim1 <- unlist(output2.Sim1)

LB <- x2.sim1[seq(1, length(x2.sim1), 2)]
UB <- x2.sim1[seq(2, length(x2.sim1), 2)]

(C.I_For_Parametric_BootStrap.Sim1 <- cbind(LB,UB))

(In_C.I_Parametric.Sim1 <-sum(C.I_For_Parametric_BootStrap.Sim1[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim1[,"UB"]>2.99999999))

(Parametric.Prop.sim1 <- mean(C.I_For_Parametric_BootStrap.Sim1[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim1[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Parametric.Sim1 <- mean(C.I_For_Parametric_BootStrap.Sim1[,"UB"]-C.I_For_Parametric_BootStrap.Sim1[,"LB"]))


#####################################################
##########' PART THREE: Re-sampling Errors ##########
#####################################################

output3.Sim1 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps = rnorm(n)
  beta0=1
  beta1=3
  y.3 = beta0+beta1*x+eps
  L3 <- lm(y.3 ~ x) 
  summary(L3)
  b03 <- L3$coeff[1]
  b13 <- L3$coeff[2]
  
  m.resid3 <- rstandard(L3, sd = 1)
  r3 <- m.resid3 - mean(m.resid3)
  
  estsErr <- replicate(m, expr={
    rstar <- sample(r3, replace=TRUE, size=n)
    ystar3 <- L3$fitted.values + rstar
    xstar3 <- x
    Lb3 <- lm(ystar3 ~ xstar3)
    s3 <- summary(Lb3)$sigma
    c(b0=Lb3$coeff[1], b1=Lb3$coeff[2], s=s3)
  })
  
  ests <- t(estsErr)
  
  slope.boot.replicates3<- ests[,2]
  
  meanCases3 <- colMeans(estsErr)
  sdCases3 <- apply(estsErr, 2, "sd")
  biasInt3 <- mean(estsErr[,1] - b03)
  biasSlope3 <- mean(estsErr[,2] - b13)
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.error.percentile <- quantile(slope.boot.replicates3,probs = p.lower,names = F)
  
  ub.error.percentile <- quantile(slope.boot.replicates3,probs = p.upper,names = F)
  
  LB <- as.vector(lb.error.percentile)
  
  output3.Sim1 <- data.frame(LB)
  
  UB <- as.vector(ub.error.percentile)
  
  output3.Sim1$UB <- UB
  
  output3.Sim1
  
})


x3.sim1 <- unlist(output3.Sim1)

LB <- x3.sim1[seq(1, length(x3.sim1), 2)]
UB <- x3.sim1[seq(2, length(x3.sim1), 2)]

(C.I_For_ReSamp.Errors_BootStrap.Sim1 <- cbind(LB,UB))

(In_C.I_Resampling.Errors.Sim1 <- sum(C.I_For_ReSamp.Errors_BootStrap.Sim1[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim1[,"UB"]>2.99999999))

(ReSamp.Error.prop.sim1 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim1[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim1[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSamp.Error.Sim1 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim1[,"UB"]-C.I_For_ReSamp.Errors_BootStrap.Sim1[,"LB"]))










#####################################################
###########" : Method 5 Wild Boot Strap #############
#####################################################


output5.Sim1 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps = rnorm(n)
  beta0=1
  beta1=3
  y.5 = beta0+beta1*x+eps
  L5 <- lm(y.5 ~ x) 
  summary(L5)
  b05 <- L5$coeff[1]
  b15 <- L5$coeff[2]
  s <- summary(L5)$sigma
  yhat <- L5$fitted.values
  
  m.resid5 <- rstandard(L5, sd = 1)
  r5 <- m.resid5 - mean(m.resid5)
  
  out.5 <- replicate(m, expr={
    v.star <- rnorm(n, 0, sd=1)
    xstar5 = x
    ystar5 <- L5$fitted.values + r5*v.star
    Lb5 <- lm(ystar5 ~ xstar5)
    s5 <- summary(Lb5)$sigma
    c(b0=Lb5$coeff[1], b1=Lb5$coeff[2], s=s5)
  })
  
  bootCases5 <- t(out.5)
  meanCases5 <- colMeans(bootCases5)
  sdCases5 <- apply(bootCases5, 2, "sd")
  
  slope.boot.replicates5<- bootCases5[,2]
  
  biasInt5 <- mean(bootCases5[,1] - b05)
  biasSlope5 <- mean(bootCases5[,2] - b15)
  rbind(estimate=c(b05, b15), bias=c(biasInt5, biasSlope5), se=sdCases5[1:2], cv=c(biasInt5, cv=biasSlope5)/sdCases5[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.wild.percentile <- quantile(slope.boot.replicates5,probs = p.lower,names = F)
  
  ub.wild.percentile <- quantile(slope.boot.replicates5,probs = p.upper,names = F)
  
  LB <- as.vector(lb.wild.percentile)
  
  output5.Sim1 <- data.frame(LB)
  
  UB <- as.vector(ub.wild.percentile)
  
  output5.Sim1$UB <- UB
  
  output5.Sim1
  
})


x5.sim1 <- unlist(output5.Sim1)

LB <- x5.sim1[seq(1, length(x5.sim1), 2)]
UB <- x5.sim1[seq(2, length(x5.sim1), 2)]

(C.I_For_Wild_BootStrap.Sim1 <- cbind(LB,UB))

(In_C.I_Wild.Sim1 <- sum(C.I_For_Wild_BootStrap.Sim1[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim1[,"UB"]>2.99999999))

(wild.prop.sim1 <- mean(C.I_For_Wild_BootStrap.Sim1[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim1[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Wild.Sim1 <- mean(C.I_For_Wild_BootStrap.Sim1[,"UB"]-C.I_For_Wild_BootStrap.Sim1[,"LB"]))



###' Simulation 1 Summary

(In_C.I_Sim1 <-rbind(In_C.I_ReSamp.Sim1, In_C.I_Parametric.Sim1, In_C.I_Resampling.Errors.Sim1, In_C.I_Wild.Sim1))

(Proportions_In_CI_Sim1 <- rbind(ReSampling.prop.sim1, Parametric.Prop.sim1, ReSamp.Error.prop.sim1, wild.prop.sim1))

(Lengths.Sim1 <- rbind(average.length.C.I.ReSampling.Cases.Sim1, average.length.C.I.Parametric.Sim1, average.length.C.I.ReSamp.Error.Sim1, average.length.C.I.Wild.Sim1))






#####################################################
#####################################################
#####################################################
##########" Simulation 2: Heavy Tail Noise ##########
#####################################################
#####################################################
#####################################################


#####################################################
##############' Model ONE: ReSampling ###############
#####################################################

output1.Sim2 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps2 = rt(n,df=3)/sqrt(3)
  beta0=1
  beta1=3
  y = beta0+beta1*x+eps2
  
  L1 <- lm(y ~ x) 
  summary(L1)
  b0 <- L1$coeff[1]
  b1 <- L1$coeff[2]
  
  m <- 2000
  out <- replicate(m, expr={
    i <- sample(1:n, replace=TRUE, size=n)
    xstar <- x[i]
    ystar <- y[i]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
  })
  
  bootCases <- t(out)
  meanCases <- colMeans(bootCases)
  sdCases <- apply(bootCases, 2, "sd")
  
  slope.boot.replicates<- bootCases[,2]
  
  biasInt <- mean(bootCases[,1] - b0)
  biasSlope <- mean(bootCases[,2] - b1)
  rbind(estimate=c(b0, b1), bias=c(biasInt, biasSlope), se=sdCases[1:2], cv=c(biasInt, cv=biasSlope)/sdCases[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.resample.percentile <- quantile(slope.boot.replicates,probs = p.lower,names = F)
  ub.resample.percentile = quantile(slope.boot.replicates,probs = p.upper,names = F)
  
  LB <- as.vector(lb.resample.percentile)
  
  output1.Sim2 <- data.frame(LB)
  
  UB <- as.vector(ub.resample.percentile)
  
  output1.Sim2$UB <- UB
  
  output1.Sim2
  
})


x1.sim2 <- unlist(output1.Sim2)

LB <- x1.sim2[seq(1, length(x1.sim2), 2)]
UB <- x1.sim2[seq(2, length(x1.sim2), 2)]

(C.I_For_ReSampling_Cases_BootStrap.Sim2 <- cbind(LB,UB))

(In_C.I_ReSamp.Sim2 <- sum(C.I_For_ReSampling_Cases_BootStrap.Sim2[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim2[,"UB"]>2.99999999))

(ReSampling.prop.sim2 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim2[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim2[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSampling.Cases.Sim2 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim2[,"UB"]-C.I_For_ReSampling_Cases_BootStrap.Sim2[,"LB"]))





#####################################################
#########' Model Two: Parametric Bootstrap ##########
#####################################################


output2.Sim2 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps2 = rt(n,df=3)/sqrt(3)
  beta0=1
  beta1=3
  y.2 = beta0+beta1*x+eps2
  L2 <- lm(y.2 ~ x)
  summary(L2)
  b0.2 <- L2$coeff[1]
  b1.2 <- L2$coeff[2]
  s <- summary(L2)$sigma
  yhat <- L2$fitted.values
  m <- 2000
  
  out.2 <- replicate(m, expr={
    eps.star.2 <- rnorm(n, 0, sd=s)
    xstar.2 <- x
    ystar.2 <- yhat+ eps.star.2
    Lb.2 <- lm(ystar.2 ~ xstar.2)
    s.2 <- summary(Lb.2)$sigma
    c(Lb.2$coeff[1], slope=Lb.2$coeff[2], s=s.2)
  })
  
  bootCases2 <- t(out.2)
  meanCases2 <- colMeans(bootCases2)
  sdCases2 <- apply(bootCases2, 2, "sd")
  
  slope.boot.replicates2<- bootCases2[,2]
  
  biasInt2 <- mean(bootCases2[,1] - b0.2)
  biasSlope2 <- mean(bootCases2[,2] - b1.2)
  rbind(estimate=c(b0.2, b1.2), bias=c(biasInt2, biasSlope2), se=sdCases2[1:2], cv=c(biasInt2, cv=biasSlope2)/sdCases2[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.parametric.percentile = quantile(slope.boot.replicates2,probs = p.lower,names = F)
  ub.parametric.percentile = quantile(slope.boot.replicates2,probs = p.upper,names = F)
  
  
  LB <- as.vector(lb.parametric.percentile)
  
  output2.Sim2 <- data.frame(LB)
  
  UB <- as.vector(ub.parametric.percentile)
  
  output2.Sim2$UB <- UB
  
  output2.Sim2
  
})


x2.sim2 <- unlist(output2.Sim2)

LB <- x2.sim2[seq(1, length(x2.sim2), 2)]
UB <- x2.sim2[seq(2, length(x2.sim2), 2)]

(C.I_For_Parametric_BootStrap.Sim2 <- cbind(LB,UB))

(In_C.I_Parametric.Sim2 <-sum(C.I_For_Parametric_BootStrap.Sim2[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim2[,"UB"]>2.99999999))

(Parametric.Prop.sim2 <- mean(C.I_For_Parametric_BootStrap.Sim2[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim2[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Parametric.Sim2 <- mean(C.I_For_Parametric_BootStrap.Sim2[,"UB"]-C.I_For_Parametric_BootStrap.Sim2[,"LB"]))



#####################################################
##########' PART THREE: Re-sampling Errors ##########
#####################################################

output3.Sim2 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps2 = rt(n,df=3)/sqrt(3)
  beta0=1
  beta1=3
  y.3 = beta0+beta1*x+eps2
  L3 <- lm(y.3 ~ x) 
  summary(L3)
  b03 <- L3$coeff[1]
  b13 <- L3$coeff[2]
  
  m.resid3 <- rstandard(L3, sd = 1)
  r3 <- m.resid3 - mean(m.resid3)
  
  estsErr <- replicate(m, expr={
    rstar <- sample(r3, replace=TRUE, size=n)
    ystar3 <- L3$fitted.values + rstar
    xstar3 <- x
    Lb3 <- lm(ystar3 ~ xstar3)
    s3 <- summary(Lb3)$sigma
    c(b0=Lb3$coeff[1], b1=Lb3$coeff[2], s=s3)
  })
  
  ests <- t(estsErr)
  
  slope.boot.replicates3<- ests[,2]
  
  meanCases3 <- colMeans(estsErr)
  sdCases3 <- apply(estsErr, 2, "sd")
  biasInt3 <- mean(estsErr[,1] - b03)
  biasSlope3 <- mean(estsErr[,2] - b13)
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.error.percentile = quantile(slope.boot.replicates3,probs = p.lower,names = F)
  ub.error.percentile = quantile(slope.boot.replicates3,probs = p.upper,names = F)
  
  
  LB <- as.vector(lb.error.percentile)
  
  output3.Sim2 <- data.frame(LB)
  
  UB <- as.vector(ub.error.percentile)
  
  output3.Sim2$UB <- UB
  
  output3.Sim2
  
})


x3.sim2 <- unlist(output3.Sim2)

LB <- x3.sim2[seq(1, length(x3.sim2), 2)]
UB <- x3.sim2[seq(2, length(x3.sim2), 2)]

(C.I_For_ReSamp.Errors_BootStrap.Sim2 <- cbind(LB,UB))

(In_C.I_Resampling.Errors.Sim2 <- sum(C.I_For_ReSamp.Errors_BootStrap.Sim2[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim2[,"UB"]>2.99999999))

(ReSamp.Error.prop.sim2 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim2[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim2[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSamp.Error.Sim2 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim2[,"UB"]-C.I_For_ReSamp.Errors_BootStrap.Sim2[,"LB"]))











#####################################################
###########" : Method 5 Wild Boot Strap #############
#####################################################


output5.Sim2 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps2 = rt(n,df=3)/sqrt(3)
  beta0=1
  beta1=3
  y.5 = beta0+beta1*x+eps2
  L5 <- lm(y.5 ~ x) 
  summary(L5)
  b05 <- L5$coeff[1]
  b15 <- L5$coeff[2]
  s <- summary(L5)$sigma
  yhat <- L5$fitted.values
  
  m.resid5 <- rstandard(L5, sd = 1)
  r5 <- m.resid5 - mean(m.resid5)
  
  out.5 <- replicate(m, expr={
    v.star <- rnorm(n, 0, sd=1)
    xstar5 = x
    ystar5 <- L5$fitted.values + r5*v.star
    Lb5 <- lm(ystar5 ~ xstar5)
    s5 <- summary(Lb5)$sigma
    c(b0=Lb5$coeff[1], b1=Lb5$coeff[2], s=s5)
  })
  
  bootCases5 <- t(out.5)
  meanCases5 <- colMeans(bootCases5)
  sdCases5 <- apply(bootCases5, 2, "sd")
  
  slope.boot.replicates5<- bootCases5[,2]
  
  biasInt5 <- mean(bootCases5[,1] - b05)
  biasSlope5 <- mean(bootCases5[,2] - b15)
  rbind(estimate=c(b05, b15), bias=c(biasInt5, biasSlope5), se=sdCases5[1:2], cv=c(biasInt5, cv=biasSlope5)/sdCases5[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  
  lb.wild.percentile = quantile(slope.boot.replicates5,probs = p.lower,names = F)
  ub.wild.percentile = quantile(slope.boot.replicates5,probs = p.upper,names = F)
  
  LB <- as.vector(lb.wild.percentile)
  
  output5.Sim2 <- data.frame(LB)
  
  UB <- as.vector(ub.wild.percentile)
  
  output5.Sim2$UB <- UB
  
  output5.Sim2
  
})


x5.sim2 <- unlist(output5.Sim2)

LB <- x5.sim2[seq(1, length(x5.sim2), 2)]
UB <- x5.sim2[seq(2, length(x5.sim2), 2)]

(C.I_For_Wild_BootStrap.Sim2 <- cbind(LB,UB))

(In_C.I_Wild.Sim2 <- sum(C.I_For_Wild_BootStrap.Sim2[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim2[,"UB"]>2.99999999))

(wild.prop.sim2 <- mean(C.I_For_Wild_BootStrap.Sim2[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim2[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Wild.Sim2 <- mean(C.I_For_Wild_BootStrap.Sim2[,"UB"]-C.I_For_Wild_BootStrap.Sim2[,"LB"]))


###' Simulation 2 Summary

(In_C.I_Sim2 <-rbind(In_C.I_ReSamp.Sim2, In_C.I_Parametric.Sim2, In_C.I_Resampling.Errors.Sim2, In_C.I_Wild.Sim2))

(Proportions_In_CI_Sim2 <- rbind(ReSampling.prop.sim2, Parametric.Prop.sim2, ReSamp.Error.prop.sim2, wild.prop.sim2))

(Lengths.Sim2 <- rbind(average.length.C.I.ReSampling.Cases.Sim2, average.length.C.I.Parametric.Sim2, average.length.C.I.ReSamp.Error.Sim2, average.length.C.I.Wild.Sim2))







#####################################################
#####################################################
#####################################################
############" Simulation 3: Skewed Noise ############
#####################################################
#####################################################
#####################################################


#####################################################
##############' Model ONE: Resampling ###############
#####################################################

output1.Sim3 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps3 = rgamma(n,shape=4,rate=2)-2
  beta0=1
  beta1=3
  y = beta0+beta1*x+eps3
  
  L1 <- lm(y ~ x) 
  summary(L1)
  b0 <- L1$coeff[1]
  b1 <- L1$coeff[2]
  
  m <- 2000
  out <- replicate(m, expr={
    i <- sample(1:n, replace=TRUE, size=n)
    xstar <- x[i]
    ystar <- y[i]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
  })
  
  bootCases <- t(out)
  meanCases <- colMeans(bootCases)
  sdCases <- apply(bootCases, 2, "sd")
  
  slope.boot.replicates<- bootCases[,2]
  
  biasInt <- mean(bootCases[,1] - b0)
  biasSlope <- mean(bootCases[,2] - b1)
  rbind(estimate=c(b0, b1), bias=c(biasInt, biasSlope), se=sdCases[1:2], cv=c(biasInt, cv=biasSlope)/sdCases[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  
  lb.resample.percentile <- quantile(slope.boot.replicates,probs = p.lower,names = F)
  ub.resample.percentile = quantile(slope.boot.replicates,probs = p.upper,names = F)
  
  LB <- as.vector(lb.resample.percentile)
  
  output1.Sim3 <- data.frame(LB)
  
  UB <- as.vector(ub.resample.percentile)
  
  output1.Sim3$UB <- UB
  
  output1.Sim3
  
})


x1.sim3 <- unlist(output1.Sim3)

LB <- x1.sim3[seq(1, length(x1.sim3), 2)]
UB <- x1.sim3[seq(2, length(x1.sim3), 2)]

(C.I_For_ReSampling_Cases_BootStrap.Sim3 <- cbind(LB,UB))

(In_C.I_ReSamp.Sim3 <- sum(C.I_For_ReSampling_Cases_BootStrap.Sim3[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim3[,"UB"]>2.99999999))

(ReSampling.prop.sim3 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim3[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim3[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSampling.Cases.Sim3 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim3[,"UB"]-C.I_For_ReSampling_Cases_BootStrap.Sim3[,"LB"]))



#####################################################
#########' Model Two: Parametric Bootstrap ##########
#####################################################

output2.Sim3 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  eps3 = rgamma(n,shape=4,rate=2)-2
  beta0=1
  beta1=3
  y.2 = beta0+beta1*x+eps3
  L2 <- lm(y.2 ~ x)
  summary(L2)
  b0.2 <- L2$coeff[1]
  b1.2 <- L2$coeff[2]
  s <- summary(L2)$sigma
  yhat <- L2$fitted.values
  m <- 2000
  
  out.2 <- replicate(m, expr={
    eps.star.2 <- rnorm(n, 0, sd=s)
    xstar.2 <- x
    ystar.2 <- yhat+ eps.star.2
    Lb.2 <- lm(ystar.2 ~ xstar.2)
    s.2 <- summary(Lb.2)$sigma
    c(Lb.2$coeff[1], slope=Lb.2$coeff[2], s=s.2)
  })
  
  bootCases2 <- t(out.2)
  meanCases2 <- colMeans(bootCases2)
  sdCases2 <- apply(bootCases2, 2, "sd")
  
  slope.boot.replicates2<- bootCases2[,2]
  
  biasInt2 <- mean(bootCases2[,1] - b0.2)
  biasSlope2 <- mean(bootCases2[,2] - b1.2)
  rbind(estimate=c(b0.2, b1.2), bias=c(biasInt2, biasSlope2), se=sdCases2[1:2], cv=c(biasInt2, cv=biasSlope2)/sdCases2[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.parametric.percentile = quantile(slope.boot.replicates2,probs = p.lower,names = F)
  ub.parametric.percentile = quantile(slope.boot.replicates2,probs = p.upper,names = F)
  
  LB <- as.vector(lb.parametric.percentile)
  
  output2.Sim3 <- data.frame(LB)
  
  UB <- as.vector(ub.parametric.percentile)
  
  output2.Sim3$UB <- UB
  
  output2.Sim3
  
})


x2.sim3 <- unlist(output2.Sim3)

LB <- x2.sim3[seq(1, length(x2.sim3), 2)]
UB <- x2.sim3[seq(2, length(x2.sim3), 2)]

(C.I_For_Parametric_BootStrap.Sim3 <- cbind(LB,UB))

(In_C.I_Parametric.Sim3 <-sum(C.I_For_Parametric_BootStrap.Sim3[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim3[,"UB"]>2.99999999))

(Parametric.Prop.sim3 <- mean(C.I_For_Parametric_BootStrap.Sim3[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim3[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Parametric.Sim3 <- mean(C.I_For_Parametric_BootStrap.Sim3[,"UB"]-C.I_For_Parametric_BootStrap.Sim3[,"LB"]))









#####################################################
##########' PART THREE: Re-sampling Errors ##########
#####################################################

output3.Sim3 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps3 = rgamma(n,shape=4,rate=2)-2
  beta0=1
  beta1=3
  y.3 = beta0+beta1*x+eps3
  L3 <- lm(y.3 ~ x) 
  summary(L3)
  b03 <- L3$coeff[1]
  b13 <- L3$coeff[2]
  
  m.resid3 <- rstandard(L3, sd = 1)
  r3 <- m.resid3 - mean(m.resid3)
  
  estsErr <- replicate(m, expr={
    rstar <- sample(r3, replace=TRUE, size=n)
    ystar3 <- L3$fitted.values + rstar
    xstar3 <- x
    Lb3 <- lm(ystar3 ~ xstar3)
    s3 <- summary(Lb3)$sigma
    c(b0=Lb3$coeff[1], b1=Lb3$coeff[2], s=s3)
  })
  
  ests <- t(estsErr)
  
  slope.boot.replicates3<- ests[,2]
  
  meanCases3 <- colMeans(estsErr)
  sdCases3 <- apply(estsErr, 2, "sd")
  biasInt3 <- mean(estsErr[,1] - b03)
  biasSlope3 <- mean(estsErr[,2] - b13)
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  
  lb.error.percentile = quantile(slope.boot.replicates3,probs = p.lower,names = F)
  ub.error.percentile = quantile(slope.boot.replicates3,probs = p.upper,names = F)
  
  
  LB <- as.vector(lb.error.percentile)
  
  output3.Sim3 <- data.frame(LB)
  
  UB <- as.vector(ub.error.percentile)
  
  output3.Sim3$UB <- UB
  
  output3.Sim3
  
})


x3.sim3 <- unlist(output3.Sim3)

LB <- x3.sim3[seq(1, length(x3.sim3), 2)]
UB <- x3.sim3[seq(2, length(x3.sim3), 2)]

(C.I_For_ReSamp.Errors_BootStrap.Sim3 <- cbind(LB,UB))

(In_C.I_Resampling.Errors.Sim3 <- sum(C.I_For_ReSamp.Errors_BootStrap.Sim3[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim3[,"UB"]>2.99999999))

(ReSamp.Error.prop.sim3 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim3[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim3[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSamp.Error.Sim3 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim3[,"UB"]-C.I_For_ReSamp.Errors_BootStrap.Sim3[,"LB"]))













#####################################################
###########" : Method 5 Wild Boot Strap #############
#####################################################

output5.Sim3 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  eps3 = rgamma(n,shape=4,rate=2)-2
  beta0=1
  beta1=3
  y.5 = beta0+beta1*x+eps3
  L5 <- lm(y.5 ~ x) 
  summary(L5)
  b05 <- L5$coeff[1]
  b15 <- L5$coeff[2]
  s <- summary(L5)$sigma
  yhat <- L5$fitted.values
  
  m.resid5 <- rstandard(L5, sd = 1)
  r5 <- m.resid5 - mean(m.resid5)
  
  out.5 <- replicate(m, expr={
    v.star <- rnorm(n, 0, sd=1)
    xstar5 = x
    ystar5 <- L5$fitted.values + r5*v.star
    Lb5 <- lm(ystar5 ~ xstar5)
    s5 <- summary(Lb5)$sigma
    c(b0=Lb5$coeff[1], b1=Lb5$coeff[2], s=s5)
  })
  
  bootCases5 <- t(out.5)
  meanCases5 <- colMeans(bootCases5)
  sdCases5 <- apply(bootCases5, 2, "sd")
  
  slope.boot.replicates5<- bootCases5[,2]
  
  biasInt5 <- mean(bootCases5[,1] - b05)
  biasSlope5 <- mean(bootCases5[,2] - b15)
  rbind(estimate=c(b05, b15), bias=c(biasInt5, biasSlope5), se=sdCases5[1:2], cv=c(biasInt5, cv=biasSlope5)/sdCases5[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.wild.percentile = quantile(slope.boot.replicates5,probs = p.lower,names = F)
  ub.wild.percentile = quantile(slope.boot.replicates5,probs = p.upper,names = F)
  
  LB <- as.vector(lb.wild.percentile)
  
  output5.Sim3 <- data.frame(LB)
  
  UB <- as.vector(ub.wild.percentile)
  
  output5.Sim3$UB <- UB
  
  output5.Sim3
  
})


x5.sim3 <- unlist(output5.Sim3)

LB <- x5.sim3[seq(1, length(x5.sim3), 2)]
UB <- x5.sim3[seq(2, length(x5.sim3), 2)]

(C.I_For_Wild_BootStrap.Sim3 <- cbind(LB,UB))

(In_C.I_Wild.Sim3 <- sum(C.I_For_Wild_BootStrap.Sim3[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim3[,"UB"]>2.99999999))

(wild.prop.sim3 <- mean(C.I_For_Wild_BootStrap.Sim3[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim3[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Wild.Sim3 <- mean(C.I_For_Wild_BootStrap.Sim3[,"UB"]-C.I_For_Wild_BootStrap.Sim3[,"LB"]))




###' Simulation 3 Summary

(In_C.I_Sim3 <-rbind(In_C.I_ReSamp.Sim3, In_C.I_Parametric.Sim3, In_C.I_Resampling.Errors.Sim3, In_C.I_Wild.Sim3))

(Proportions_In_CI_Sim3 <- rbind(ReSampling.prop.sim3, Parametric.Prop.sim3, ReSamp.Error.prop.sim3, wild.prop.sim3))

(Lengths.Sim3 <- rbind(average.length.C.I.ReSampling.Cases.Sim3, average.length.C.I.Parametric.Sim3, average.length.C.I.ReSamp.Error.Sim3, average.length.C.I.Wild.Sim3))




#####################################################
#####################################################
#####################################################
######" Simulation 4: Heteroscedastic Noise #########
#####################################################
#####################################################
#####################################################


#####################################################
##############' Model ONE: Resampling ###############
#####################################################


###' Lower Bound

output1.Sim4 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  esp = rnorm(n)*x^2/sqrt(mean(x^4))
  esp4 = esp-mean(esp)
  beta0=1
  beta1=3
  y = beta0+beta1*x+esp4
  
  L1 <- lm(y ~ x) 
  summary(L1)
  b0 <- L1$coeff[1]
  b1 <- L1$coeff[2]
  
  m <- 2000
  out <- replicate(m, expr={
    i <- sample(1:n, replace=TRUE, size=n)
    xstar <- x[i]
    ystar <- y[i]
    Lb <- lm(ystar ~ xstar)
    s <- summary(Lb)$sigma
    c(Lb$coeff[1], slope=Lb$coeff[2], s=s)
  })
  
  bootCases <- t(out)
  meanCases <- colMeans(bootCases)
  sdCases <- apply(bootCases, 2, "sd")
  
  slope.boot.replicates<- bootCases[,2]
  
  biasInt <- mean(bootCases[,1] - b0)
  biasSlope <- mean(bootCases[,2] - b1)
  rbind(estimate=c(b0, b1), bias=c(biasInt, biasSlope), se=sdCases[1:2], cv=c(biasInt, cv=biasSlope)/sdCases[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.resample.percentile <- quantile(slope.boot.replicates,probs = p.lower,names = F)
  
  ub.resample.percentile = quantile(slope.boot.replicates,probs = p.upper,names = F)
  
  LB <- as.vector(lb.resample.percentile)
  
  output1.Sim4 <- data.frame(LB)
  
  UB <- as.vector(ub.resample.percentile)
  
  output1.Sim4$UB <- UB
  
  output1.Sim4
  
})


x1.sim4 <- unlist(output1.Sim4)

LB <- x1.sim4[seq(1, length(x1.sim4), 2)]
UB <- x1.sim4[seq(2, length(x1.sim4), 2)]

(C.I_For_ReSampling_Cases_BootStrap.Sim4 <- cbind(LB,UB))

(In_C.I_ReSamp.Sim4 <- sum(C.I_For_ReSampling_Cases_BootStrap.Sim4[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim4[,"UB"]>2.99999999))

(ReSampling.prop.sim4 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim4[,"LB"]<3.00000001&C.I_For_ReSampling_Cases_BootStrap.Sim4[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSampling.Cases.Sim4 <- mean(C.I_For_ReSampling_Cases_BootStrap.Sim4[,"UB"]-C.I_For_ReSampling_Cases_BootStrap.Sim4[,"LB"]))









#####################################################
#########' Model Two: Parametric Bootstrap ##########
#####################################################

output2.Sim4 = replicate(w, expr={
  
  n = 2000
  x = runif(n, -1, 1)
  esp = rnorm(n)*x^2/sqrt(mean(x^4))
  esp4 = esp-mean(esp)
  beta0=1
  beta1=3
  y.2 = beta0+beta1*x+esp4
  L2 <- lm(y.2 ~ x)
  summary(L2)
  b0.2 <- L2$coeff[1]
  b1.2 <- L2$coeff[2]
  s <- summary(L2)$sigma
  yhat <- L2$fitted.values
  m <- 2000
  
  out.2 <- replicate(m, expr={
    eps.star.2 <- rnorm(n, 0, sd=s)
    xstar.2 <- x
    ystar.2 <- yhat+ eps.star.2
    Lb.2 <- lm(ystar.2 ~ xstar.2)
    s.2 <- summary(Lb.2)$sigma
    c(Lb.2$coeff[1], slope=Lb.2$coeff[2], s=s.2)
  })
  
  bootCases2 <- t(out.2)
  meanCases2 <- colMeans(bootCases2)
  sdCases2 <- apply(bootCases2, 2, "sd")
  
  slope.boot.replicates2<- bootCases2[,2]
  
  biasInt2 <- mean(bootCases2[,1] - b0.2)
  biasSlope2 <- mean(bootCases2[,2] - b1.2)
  rbind(estimate=c(b0.2, b1.2), bias=c(biasInt2, biasSlope2), se=sdCases2[1:2], cv=c(biasInt2, cv=biasSlope2)/sdCases2[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.parametric.percentile = quantile(slope.boot.replicates2,probs = p.lower,names = F)
  
  ub.parametric.percentile = quantile(slope.boot.replicates2,probs = p.upper,names = F)
  
  LB <- as.vector(lb.parametric.percentile)
  
  output2.Sim4 <- data.frame(LB)
  
  UB <- as.vector(ub.parametric.percentile)
  
  output2.Sim4$UB <- UB
  
  output2.Sim4
  
})


x2.sim4 <- unlist(output2.Sim4)

LB <- x2.sim4[seq(1, length(x2.sim4), 2)]
UB <- x2.sim4[seq(2, length(x2.sim4), 2)]

(C.I_For_Parametric_BootStrap.Sim4 <- cbind(LB,UB))

(In_C.I_Parametric.Sim4 <-sum(C.I_For_Parametric_BootStrap.Sim4[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim4[,"UB"]>2.99999999))

(Parametric.Prop.sim4 <- mean(C.I_For_Parametric_BootStrap.Sim4[,"LB"]<3.00000001&C.I_For_Parametric_BootStrap.Sim4[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Parametric.Sim4 <- mean(C.I_For_Parametric_BootStrap.Sim4[,"UB"]-C.I_For_Parametric_BootStrap.Sim4[,"LB"]))







#####################################################
##########' PART THREE: Re-sampling Errors ##########
#####################################################

output3.Sim4 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  esp = rnorm(n)*x^2/sqrt(mean(x^4))
  esp4 = esp-mean(esp)
  beta0=1
  beta1=3
  y.3 = beta0+beta1*x+esp4
  L3 <- lm(y.3 ~ x) 
  summary(L3)
  b03 <- L3$coeff[1]
  b13 <- L3$coeff[2]
  
  m.resid3 <- rstandard(L3, sd = 1)
  r3 <- m.resid3 - mean(m.resid3)
  
  estsErr <- replicate(m, expr={
    rstar <- sample(r3, replace=TRUE, size=n)
    ystar3 <- L3$fitted.values + rstar
    xstar3 <- x
    Lb3 <- lm(ystar3 ~ xstar3)
    s3 <- summary(Lb3)$sigma
    c(b0=Lb3$coeff[1], b1=Lb3$coeff[2], s=s3)
  })
  
  ests <- t(estsErr)
  
  slope.boot.replicates3<- ests[,2]
  
  meanCases3 <- colMeans(estsErr)
  sdCases3 <- apply(estsErr, 2, "sd")
  biasInt3 <- mean(estsErr[,1] - b03)
  biasSlope3 <- mean(estsErr[,2] - b13)
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.error.percentile = quantile(slope.boot.replicates3,probs = p.lower,names = F)
  
  ub.error.percentile = quantile(slope.boot.replicates3,probs = p.upper,names = F)
  
  LB <- as.vector(lb.error.percentile)
  
  output3.Sim4 <- data.frame(LB)
  
  UB <- as.vector(ub.error.percentile)
  
  output3.Sim4$UB <- UB
  
  output3.Sim4
  
})


x3.sim4 <- unlist(output3.Sim4)

LB <- x3.sim4[seq(1, length(x3.sim4), 2)]
UB <- x3.sim4[seq(2, length(x3.sim4), 2)]

(C.I_For_ReSamp.Errors_BootStrap.Sim4 <- cbind(LB,UB))

(In_C.I_Resampling.Errors.Sim4 <- sum(C.I_For_ReSamp.Errors_BootStrap.Sim4[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim4[,"UB"]>2.99999999))

(ReSamp.Error.prop.sim4 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim4[,"LB"]<3.0001&C.I_For_ReSamp.Errors_BootStrap.Sim4[,"UB"]>2.99999999))

###' Length

(average.length.C.I.ReSamp.Error.Sim4 <- mean(C.I_For_ReSamp.Errors_BootStrap.Sim4[,"UB"]-C.I_For_ReSamp.Errors_BootStrap.Sim4[,"LB"]))



#####################################################
###########" : Method 4 Wild Boot Strap #############
#####################################################


output5.Sim4 = replicate(w, expr={
  
  m <- 2000
  n <- 2000
  x = runif(n, -1, 1)
  esp = rnorm(n)*x^2/sqrt(mean(x^4))
  esp4 = esp-mean(esp)
  beta0=1
  beta1=3
  y.5 = beta0+beta1*x+esp4
  L5 <- lm(y.5 ~ x) 
  summary(L5)
  b05 <- L5$coeff[1]
  b15 <- L5$coeff[2]
  s <- summary(L5)$sigma
  yhat <- L5$fitted.values
  
  m.resid5 <- rstandard(L5, sd = 1)
  r5 <- m.resid5 - mean(m.resid5)
  
  out.5 <- replicate(m, expr={
    v.star <- rnorm(n, 0, sd=1)
    xstar5 = x
    ystar5 <- L5$fitted.values + r5*v.star
    Lb5 <- lm(ystar5 ~ xstar5)
    s5 <- summary(Lb5)$sigma
    c(b0=Lb5$coeff[1], b1=Lb5$coeff[2], s=s5)
  })
  
  bootCases5 <- t(out.5)
  meanCases5 <- colMeans(bootCases5)
  sdCases5 <- apply(bootCases5, 2, "sd")
  
  slope.boot.replicates5<- bootCases5[,2]
  
  biasInt5 <- mean(bootCases5[,1] - b05)
  biasSlope5 <- mean(bootCases5[,2] - b15)
  rbind(estimate=c(b05, b15), bias=c(biasInt5, biasSlope5), se=sdCases5[1:2], cv=c(biasInt5, cv=biasSlope5)/sdCases5[1:2])
  
  p.upper = 1-0.1/2
  p.lower=0.1/2
  
  lb.wild.percentile = quantile(slope.boot.replicates5,probs = p.lower,names = F)
  ub.wild.percentile = quantile(slope.boot.replicates5,probs = p.upper,names = F)
  
  LB <- as.vector(lb.wild.percentile)
  
  output5.Sim4 <- data.frame(LB)
  
  UB <- as.vector(ub.wild.percentile)
  
  output5.Sim4$UB <- UB
  
  output5.Sim4
  
})


x5.sim4 <- unlist(output5.Sim4)

LB <- x5.sim4[seq(1, length(x5.sim4), 2)]
UB <- x5.sim4[seq(2, length(x5.sim4), 2)]

(C.I_For_Wild_BootStrap.Sim4 <- cbind(LB,UB))

(In_C.I_Wild.Sim4 <- sum(C.I_For_Wild_BootStrap.Sim4[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim4[,"UB"]>2.99999999))

(wild.prop.sim4 <- mean(C.I_For_Wild_BootStrap.Sim4[,"LB"]<3.0000001&C.I_For_Wild_BootStrap.Sim4[,"UB"]>2.99999999))

###' Length

(average.length.C.I.Wild.Sim4 <- mean(C.I_For_Wild_BootStrap.Sim4[,"UB"]-C.I_For_Wild_BootStrap.Sim4[,"LB"]))





###' Simulation 4 Summary

(In_C.I_Sim4 <-rbind(In_C.I_ReSamp.Sim4, In_C.I_Parametric.Sim4, In_C.I_Resampling.Errors.Sim4, In_C.I_Wild.Sim4))

(Proportions_In_CI_Sim4 <- rbind(ReSampling.prop.sim4, Parametric.Prop.sim4, ReSamp.Error.prop.sim4, wild.prop.sim4))

(Lengths.Sim4 <- rbind(average.length.C.I.ReSampling.Cases.Sim4, average.length.C.I.Parametric.Sim4, average.length.C.I.ReSamp.Error.Sim4, average.length.C.I.Wild.Sim4))

?rnorm()


###' Overall Summary Across Simulations

cbind(Proportions_In_CI_Sim1, Proportions_In_CI_Sim2, Proportions_In_CI_Sim3, Proportions_In_CI_Sim4)

cbind(Lengths.Sim1, Lengths.Sim2, Lengths.Sim3, Lengths.Sim4)

2.929292-3.112795