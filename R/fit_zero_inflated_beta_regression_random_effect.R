
cal_zibeta_loglik = function(para,
                           subject.n,time.n,
                           X.aug,Z.aug,Y,
                           X.test.coeff.index,
                           Z.test.coeff.index,
                           prod.mat,
                           gh.weights,gh.nodes,
                           quad.n=quad.n){
  #### s1, s2, v will always be estimated
  logistic.para <- para[1:(sum(!X.test.coeff.index)+1)] ## s1, alpha
  beta.para <- para[-(1:(sum(!X.test.coeff.index)+1))] ## s2,v, beta
  logistic.logLike <- cal_logistic_loglik(para=logistic.para,
                                          X.aug=X.aug, Y=Y,
                                          subject.n=subject.n,
                                          time.n=time.n,
                                          prod.mat=prod.mat,
                                          X.test.coeff.index=X.test.coeff.index,
                                          gh.weights=gh.weights,gh.nodes=gh.nodes,
                                          quad.n=quad.n)
  beta.logLike <- cal_beta_loglik(para=beta.para,
                                  Z.aug=Z.aug,Y=Y,
                                  subject.n=subject.n,
                                  time.n=time.n,
                                  Z.test.coeff.index=Z.test.coeff.index, 
                                  prod.mat=prod.mat,
                                  gh.weights=gh.weights,gh.nodes=gh.nodes,
                                  quad.n=quad.n)
  return(logistic.logLike+beta.logLike)
}



#######################################
fit_zero_inflated_beta_random_effect = function(X=X,Z=Z,Y=Y,
                                  subject.ind=subject.ind,time.ind=time.ind,
                                  quad.n=30,verbose=FALSE){
  ######
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  X.aug <- cbind(intersept = 1, X)
  Z.aug <- cbind(intersept = 1, Z)
  #############
  subject.n <- length(unique(subject.ind))
  time.n    <- length(unique(time.ind))
  prod.mat <- matrix(rep(c(rep(1,time.n),rep(0,subject.n*time.n)),subject.n)[1:(subject.n^2*time.n)],byrow=TRUE,nrow=subject.n,ncol=subject.n*time.n)
  #############
  #### generate quad points
  gherm <- generate_gaussian_quad_points(quad.n = quad.n)
  gh.weights <- matrix(rep(gherm$weights,subject.n),nrow = subject.n,byrow = TRUE)
  gh.nodes <- matrix(rep(gherm$nodes,subject.n * time.n),
                     nrow = subject.n * time.n,byrow = TRUE)
  #############
  ###### estimate and test each parameter
  logistic.fit <- fit_logistic_random_effect(X=X,Y=Y,
                            subject.ind=subject.ind,time.ind=time.ind,
                            quad.n=quad.n,verbose=verbose)
  beta.fit <- fit_beta_random_effect(Z=Z,Y=Y,subject.ind=subject.ind,time.ind=time.ind,
                         quad.n=quad.n,verbose=verbose)
  ##################################
  ##### test the overall difference
  ##### H1: estimate all parameters
  X.test.coeff.index <- rep(FALSE,ncol(X.aug))
  Z.test.coeff.index <- rep(FALSE,ncol(Z.aug))
  opt.H1 <- nlminb(start= c(c(1,rep(0,sum(!X.test.coeff.index))), ## s1,alpha
                          c(1,2,rep(0,sum(!Z.test.coeff.index)))), ## s2,v,beta
                   objective=cal_zibeta_loglik,
                   lower = c(c(0.00001,rep(-Inf,sum(!X.test.coeff.index))),
                           c(0.00001,0.00001,rep(-Inf,sum(!Z.test.coeff.index)))),
                   upper = Inf,
                   X.test.coeff.index = X.test.coeff.index,
                   Z.test.coeff.index = Z.test.coeff.index,
                   Y=Y,X.aug=X.aug,Z.aug=Z.aug,time.n=time.n,subject.n=subject.n,
                   prod.mat=prod.mat,
                   gh.weights=gh.weights,gh.nodes=gh.nodes,
                   quad.n=quad.n,
                   control=list(trace=ifelse(verbose,2,0))
  )
  ###### H0:set all regression coefficients to zero, excluding intercept
  X.test.coeff.index <- c(FALSE,rep(TRUE,ncol(X)))
  Z.test.coeff.index <- c(FALSE,rep(TRUE,ncol(Z)))
  opt.H0 <- nlminb(start= c(c(1,rep(0,sum(!X.test.coeff.index))), ## s1,alpha
                            c(1,2,rep(0,sum(!Z.test.coeff.index)))), ## s2,v,beta
                   objective=cal_zibeta_loglik,
                   lower = c(c(0.00001,rep(-Inf,sum(!X.test.coeff.index))),
                             c(0.00001,0.00001,rep(-Inf,sum(!Z.test.coeff.index)))),
                   upper = Inf,
                   X.test.coeff.index = X.test.coeff.index,
                   Z.test.coeff.index = Z.test.coeff.index,
                   Y=Y,X.aug=X.aug,Z.aug=Z.aug,time.n=time.n,subject.n=subject.n,
                   prod.mat=prod.mat,
                   gh.weights=gh.weights,gh.nodes=gh.nodes,
                   quad.n=quad.n,
                   control=list(trace=ifelse(verbose,2,0))
  )
  likelihodd.ratio <- -2*(-opt.H0$objective-(-opt.H1$objective))
  overall.LRT.p <- 1-pchisq(likelihodd.ratio,df=sum(X.test.coeff.index)+sum(Z.test.coeff.index))

  return(list(logistic.est.table=logistic.fit$est.table,
              logistic.s1.est=logistic.fit$s1.est,
              beta.est.table=beta.fit$est.table,
              beta.s2.est=beta.fit$s2.est,
              beta.v.est=beta.fit$v.est,
              overall.test.p=overall.LRT.p))
}
