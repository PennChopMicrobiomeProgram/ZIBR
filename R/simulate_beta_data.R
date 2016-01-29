simulate_beta_random_effect_data = function(subject.n=50, time.n=5, v=2,
                              beta=as.matrix(c(-0.5,-0.5,0.5)),
                              Z=NA, s2=1, sim.seed=100){
  ##########################
  ##########################
  #### Beta regression
  if (length(NA)==1 & any(is.na(Z))){
    set.seed(sim.seed*5000+10)
    Z <- as.matrix(data.frame(log.Time=as.matrix(log(rep(seq(1,time.n),subject.n))),
            Treatment=as.matrix(c(rep(0,subject.n*time.n/2),rep(1,subject.n*time.n/2)))))
  }
  if(is.null(colnames(Z))){
    Z <- as.matrix(Z)
    colnames(Z) <- paste('var',seq(1,ncol(Z)),sep='')
  }
  
  set.seed(sim.seed*5000+1)
  c <- as.matrix(rnorm(subject.n,mean=0,sd=s2))
  c.rep <- as.matrix(as.vector(matrix(c,nrow=time.n,ncol=length(c),byrow=TRUE)))
  #####
  subject.ind <- as.vector(matrix(paste('Subject_',seq(1,subject.n),sep=''),nrow=time.n,ncol=subject.n,byrow=TRUE))
  time.ind  <- rep(seq(1,time.n),subject.n)
  ######
  Z.aug <- cbind(intersept=1,Z)
  #browser()
  logit.u  <- Z.aug %*% beta + c.rep
  ######
  u <- 1 / (1 + exp(-logit.u))
  ## v is the phi

  ######
  set.seed(sim.seed*5000+4)
  Y <- rbeta(subject.n*time.n, shape1 = u*v, shape2=(1-u)*v)
  if(any(Y>1-10^(-6))){Y[Y>1-10^(-6)] <- 1-10^(-6)}
  
  #### For test purpose
  #### betareg can not fit random effect model
  #### so set the s2 to a small value (small random effect)
  #library(betareg)
  #tdata <- data.frame(Y=Y,Z,SID=subject.ind)
  #gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='ML')
  #summary(gy)
  #gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BC')
  #summary(gy)
  #gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BR')
  #summary(gy)
  #
  #
  #fit_beta_random_effect(Z=Z,Y=Y,subject.ind=subject.ind,time.ind=time.ind,
  #                      quad.n=30,verbose=FALSE)
  
  return(list(Y=Y,Z=Z,c=c,u=u,v=v,beta=beta,s2=s2,
              subject.ind=subject.ind,time.ind=time.ind))
}
