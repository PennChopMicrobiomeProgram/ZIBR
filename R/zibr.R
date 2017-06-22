#' Fit zero-inflated beta regression with random effects
#'
#' @param logistic.cov the covariates in logistic component
#' @param beta.cov the covariates in beta component
#' @param Y the response variable in the regression model
#' @param subject.ind the variable for subject IDs
#' @param time.ind the variable for time points
#' @param component.wise.test whether to perform component wise test. If true, ZIBR will calculate pvalues for logistic and beta component respectively.
#' @param quad.n Gaussian quadrature points
#' @param verbose print the fitting process
#' @return logistic.est.table the estimated coefficients for logistic component.
#' @return logistic.s1.est the estimated standard deviation for the random effect in the logistic component.
#' @return beta.est.table the estimated coefficients for logistic component.
#' @return beta.s2.est the estimated standard deviation for the random effect in the beta component.
#' @return beta.v.est the estiamted dispersion parameter in the beta component.
#' @return loglikelihood  the log likelihood of fitting zibr model on the data.
#' @return joint.p  the pvalues for jointly testing each covariate in both logistic and beta component.
#' @export
#' @examples
#' \dontrun{
#' ## simulate some data
#' sim <- simulate_zero_inflated_beta_random_effect_data(
#' subject.n=100,time.n=5,
#' X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#' Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
#' alpha = as.matrix(c(-0.5,1)),
#' beta = as.matrix(c(-0.5,0.5)),
#' s1 = 1,s2 = 0.8,
#' v = 5,
#' sim.seed=100)
#' ## run zibr on the simulated data
#' zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$Z, Y = sim$Y,
#'                  subject.ind = sim$subject.ind,time.ind = sim$time.ind)
#' zibr.fit
#' }




zibr = function(logistic.cov=logistic.cov,
                 beta.cov=beta.cov,
                 Y=Y,
                 subject.ind=subject.ind,
                 time.ind=time.ind,
                 component.wise.test=TRUE,
                 quad.n=30,verbose=FALSE){
  Y <- as.matrix(Y)
  #### check if Y is in [0,1)
  if(min(Y)<0 | max(Y)>=1){
    stop("The response variable should be in [0,1)")
  }
  #### check the dimentions of X,Z,Y
  if(nrow(logistic.cov) != nrow(Y) | nrow(beta.cov) != nrow(Y) ){
    stop("The dimensions of covariates and repsonse variable are not correct")
  }
  #### check how many zeros in Y
  if(sum(Y>0)/length(Y)>0.9){warning("Too few zeros in the abundance data. The logistic component may be not accurate.")}
  if(sum(Y>0)/length(Y)<0.1){warning("Too many zeros in the abundance data. The beta component may be not accurate.")}
  #### if the colnames are the same, jointly test the two component
  if (identical(colnames(logistic.cov),colnames(beta.cov))){
    joint.test <- TRUE
  } 
  else{
    joint.test <- FALSE
  }
  #### check if time.ind are the same for each subject.ind
  fit = fit_zero_inflated_beta_random_effect(X=logistic.cov,Z=beta.cov,Y=Y,
 subject.ind=subject.ind,time.ind=time.ind,joint.test=joint.test)

  return(list(logistic.est.table=fit$logistic.est.table,
              logistic.s1.est=fit$logistic.s1.est,
              beta.est.table=fit$beta.est.table,
              beta.s2.est=fit$beta.s2.est,
              beta.v.est=fit$beta.v.est,
              loglikelihood = fit$loglikelihood,
              joint.p=fit$joint.p))

}
