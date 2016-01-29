
simulate_zero_inflated_beta_random_effect_data <- function(
                               subject.n=50, time.n=5, v=2,
                               alpha=as.matrix(c(0,0.5,-1)),
                               beta=as.matrix(c(-0.5,-0.5,0.5)),
                               X=NA,Z=NA,
                               s1=0.2, s2=0.2, sim.seed=100){

  ######
  if (length(NA)==1 & any(is.na(X))){
    set.seed(sim.seed*5000+10)
    X <- as.matrix(data.frame(log.Time=as.matrix(log(rep(seq(1,time.n),subject.n))),
        Treatment=as.matrix(c(rep(0,subject.n*time.n/2),rep(1,subject.n*time.n/2)))))
    #X <- as.matrix(runif(subject.n*time.n,-1,1))
    #X  <- as.matrix(c(rep(0,N*time.n/2),rep(1,N*time.n/2)))
  }
  if(is.null(colnames(X))){
    X <- as.matrix(X)
    colnames(X) <- paste('var',seq(1,ncol(X)),sep='')
  }
  if (length(NA)==1 & any(is.na(Z))){
    Z <- X
  }
  
  
  set.seed(sim.seed*5000+1)
  b <- as.matrix(rnorm(subject.n,mean=0,sd=s1))
  b.rep <- as.matrix(as.vector(matrix(b,nrow=time.n,ncol=length(b),byrow=TRUE)))
  set.seed(sim.seed*5000+2)
  c <- as.matrix(rnorm(subject.n,mean=0,sd=s2))
  c.rep <- as.matrix(as.vector(matrix(c,nrow=time.n,ncol=length(c),byrow=TRUE)))
  #####
  subject.ind <- as.vector(matrix(paste('Subject_',seq(1,subject.n),sep=''),nrow=time.n,ncol=subject.n,byrow=TRUE))
  time.ind  <- rep(seq(1,time.n),subject.n)
  ######
  X.aug <- cbind(interespt=1,X)
  Z.aug <- cbind(intersept=1,Z)
  #browser()
  logit.p  <- X.aug %*% alpha + b.rep
  logit.u  <- Z.aug %*% beta + c.rep
  #print(Z.aug %*% beta)
  ##### beta can not be too big, otherwise u will be close to 1

  ######
  p  <- 1 / (1 + exp(-logit.p))
  u <- 1 / (1 + exp(-logit.u))
  #set.seed(sim.seed+2)
  #v <- runif(1,0,5) ## v is the phi
  if (is.na(v)){v <- 2}


  ######
  Y = rep(NA,subject.n*time.n)
  set.seed(sim.seed*5000+3)
  ind.mix <- rbinom(length(Y), 1, p)
  for (i in 1:length(Y)){
    if(ind.mix[i] == 0){Y[i]=0}
    if(ind.mix[i] == 1){
      #rbeta(n, shape1, shape2, ncp = 0)
      set.seed(sim.seed*5000+i)
      Y[i] = rbeta(1, shape1 = u[i]*v, shape2=(1-u[i])*v)
      if(Y[i]>1-10^(-6)){Y[i]=1-10^(-6)}
      #while(round(Y[i],6)==1){
      #  Y[i] = rbeta(1, shape1 = u[i]*v, shape2=(1-u[i])*v)
      #}
    }
  }

  return(list(Y=Y,X=X,Z=Z,b=b,c=c,u=u,v=v,alpha=alpha,beta=beta,s1=s1,s2=s2,
              subject.ind=subject.ind,time.ind=time.ind))
}
