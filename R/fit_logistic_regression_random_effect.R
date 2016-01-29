
cal_logistic_loglik = function(para,X.aug,Y,subject.n,time.n,
                               prod.mat,
                               X.test.coeff.index,
                               gh.weights,gh.nodes,quad.n){
  ######################################
  #### given the values of random effect b
  #### calculate the probablity p
  ## X.aug is (subject.n*time.n) X (number of covariates + 1)
  ## Y and p are (subject.n*time.n) X 1
  ## b is (subject.n*time.n) X (number of quadurate points)
  s1 <- para[1]
  #### X.test.coeff.index == TRUE means we want to test this parameter
  #### Therefore we want to keep it as 0 for LRT
  alpha <- rep(NA,length(X.test.coeff.index))
  alpha[X.test.coeff.index]  <- 0 ## initialized as 0
  alpha[!X.test.coeff.index] <- para[-1][1:sum(!X.test.coeff.index)]
  alpha <- as.matrix(alpha)
  #### p can not be 0 or 1. Due to round error, logL will be NA
  #### manually set p=0 or p=1 to some other values
  p  <- 1 / (1 + exp(-(X.aug%*%alpha[,rep(1,quad.n)] + gh.nodes*s1*sqrt(2))) )
  p[p<10^(-7)] <- 10^(-7)
  p[p>(1-10^(-7))] <- (1-10^(-7))
  ### because p^[I(Y>0] * (1-p)^[I[Y=0]], replace p[Y==0] with 1-p[Y==0]
  p[Y == 0,] <- 1 - p[Y == 0,]
  #### exp(log(A*B)) = exp(logA+logB)=A*B
  logL <- sum(log(rowSums(gh.weights / sqrt(pi) * exp(prod.mat %*% log(p)))))
  return(-logL)
}

### b <- array(c(1:250000, 1:250000),c(5000,5000,2))
### system.time(rs4 <- colSums(aperm(b, c(2,1,3))))

#######################################
fit_logistic_random_effect = function(X=X,Y=Y,
                                subject.ind=subject.ind,time.ind=time.ind,
                                quad.n=30,verbose=FALSE){

  ######
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  X.aug <- cbind(intersept = 1, X)
  ######
  est.table <- matrix(NA,ncol=2,nrow=ncol(X.aug),
                      dimnames=list(colnames(X.aug),c('Estimate','Pvalue')))
  #############
  subject.n <- length(unique(subject.ind))
  time.n   <- length(unique(time.ind))
  prod.mat <- matrix(rep(c(rep(1,time.n),rep(0,subject.n*time.n)),subject.n)[1:(subject.n^2*time.n)],byrow=TRUE,nrow=subject.n,ncol=subject.n*time.n)
  #############
  #### generate quad points
  gherm <- generate_gaussian_quad_points(quad.n = quad.n)
  gh.weights <- matrix(rep(gherm$weights,subject.n),nrow = subject.n,byrow = TRUE)
  gh.nodes <- matrix(rep(gherm$nodes,subject.n * time.n),
                     nrow = subject.n * time.n,byrow = TRUE)
  #############
  #### re-order X,Y so that values belongs to the same subject are together
  #### need the values in this format for loglikelihood calculation
  gind <- sort(subject.ind,index.return=TRUE)$ix
  Y <- Y[gind,]
  X <- X[gind,]
  X.aug <- X.aug[gind,]
  ##### H1: estimate all parameters
  X.test.coeff.index <- rep(FALSE,ncol(X.aug))
  #browser()
  opt.H1 <- nlminb(start=c(1,rep(0,sum(!X.test.coeff.index))), ## s1,alpha
                   objective=cal_logistic_loglik,
                   lower = c(0.00001,
                             rep(-Inf,sum(!X.test.coeff.index))),
                   upper = Inf,
                   X.test.coeff.index = X.test.coeff.index,
                   Y=Y,X.aug=X.aug,time.n=time.n,subject.n=subject.n,
                   gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                   prod.mat=prod.mat,
                   control=list(trace=ifelse(verbose,2,0))
  )
  #### save the estimatd results
  s1.est <- opt.H1$par[1]
  alpha.est <- opt.H1$par[-1]
  est.table[,'Estimate'] <- alpha.est
  ########################
  ####### H0
  for (test.i in 1:ncol(X.aug)){
    X.test.coeff.index <- rep(FALSE,ncol(X.aug))
    X.test.coeff.index[test.i] <- TRUE
    opt.H0 <- nlminb(start=c(1,rep(0,sum(!X.test.coeff.index))), ## s1,alpha
                          objective=cal_logistic_loglik,
                          lower = c(0.00001,
                                    rep(-Inf,sum(!X.test.coeff.index))),
                          upper = c(Inf),
                          X.test.coeff.index = X.test.coeff.index,
                          Y=Y,X.aug=X.aug,time.n=time.n,subject.n=subject.n,
                          gh.weights=gh.weights,gh.nodes=gh.nodes,quad.n=quad.n,
                          prod.mat=prod.mat,
                          control=list(trace=ifelse(verbose,2,0))
    )


    likelihodd.ratio <- -2*(-opt.H0$objective-(-opt.H1$objective))
    LRT.p <- 1-pchisq(likelihodd.ratio,df=1)
    est.table[test.i,'Pvalue'] <- LRT.p
  }

  return(list(est.table=est.table,s1.est=s1.est))
}
