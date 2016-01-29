#' Fit zero-inflated beta regression with random effects
#'
#' @param logistic.cov the covariates in logistic component
#' @param beta.cov the covariates in beta component
#' @param Y the response variable in the regression model
#' @param subject.ind the variable for subjects
#' @param time.ind the variable for time points
#' @param quad.n Gaussian quadrature points
#' @param verbose print the fitting process
#' @return X,Y,Z
#' @export
#' @examples
#' \dontrun{
#' zibre(logistic.cov=logistic.cov,beta.cov=beta.cov,Y=Y,subject.ind=subject.ind,time.ind=time.ind)
#' }




zibre = function(logistic.cov=logistic.cov,
                 beta.cov=beta.cov,
                 Y=Y,
                 subject.ind=subject.ind,
                 time.ind=time.ind,
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
  #### check if time.ind are the same for each subject.ind
  fit_zero_inflated_beta_random_effect(X=logistic.cov,Z=beta.cov,Y=Y,
 subject.ind=subject.ind,time.ind=time.ind)

}
