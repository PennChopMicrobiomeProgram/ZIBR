#' Fit zero-inflated beta regression with random effects
#'
#' @param logistic_cov the covariates in logistic component
#' @param beta_cov the covariates in beta component
#' @param Y the response variable in the regression model
#' @param subject_ind the variable for subject IDs
#' @param time_ind the variable for time points
#' @param component_wise_test whether to perform component wise test.
#'   If true, ZIBR will calculate p-values for logistic and beta component respectively.
#' @param quad_n Gaussian quadrature points
#' @param verbose print the fitting process
#' @return a named list
#' \itemize{
#'   \item logistic_est_table - the estimated coefficients for logistic component.
#'   \item logistic_s1_est - the estimated standard deviation for the random effect in the logistic component.
#'   \item beta_est_table - the estimated coefficients for logistic component.
#'   \item beta_s2_est - the estimated standard deviation for the random effect in the beta component.
#'   \item beta_v_est - the estimated dispersion parameter in the beta component.
#'   \item loglikelihood - the log likelihood of fitting ZIBR model on the data.
#'   \item joint_p - the p-values for jointly testing each covariate in both logistic and beta component.
#' }
#'
#' @export
#' @examples
#' ## simulate some data
#' sim <- simulate_zero_inflated_beta_random_effect_data(
#'   subject_n = 100, time_n = 5,
#'   X = as.matrix(c(rep(0, 50 * 5), rep(1, 50 * 5))),
#'   Z = as.matrix(c(rep(0, 50 * 5), rep(1, 50 * 5))),
#'   alpha = as.matrix(c(-0.5, 1)),
#'   beta = as.matrix(c(-0.5, 0.5)),
#'   s1 = 1, s2 = 0.8,
#'   v = 5,
#'   sim_seed = 100
#' )
#'
#' ## run zibr on the simulated data
#' zibr_fit <- zibr(
#'   logistic_cov = sim$X, beta_cov = sim$Z, Y = sim$Y,
#'   subject_ind = sim$subject_ind, time_ind = sim$time_ind
#' )
#'
#' zibr_fit
zibr <- function(logistic_cov,
                 beta_cov,
                 Y,
                 subject_ind,
                 time_ind,
                 component_wise_test = TRUE,
                 quad_n = 30,
                 verbose = FALSE) {
  Y <- as.matrix(Y)
  #### check if Y is in [0,1)
  if (min(Y) < 0 || max(Y) >= 1) {
    stop("The response variable should be in [0,1)")
  }
  #### check the dimentions of X,Z,Y
  if (nrow(logistic_cov) != nrow(Y) || nrow(beta_cov) != nrow(Y)) {
    stop("The dimensions of covariates and repsonse variable are not correct")
  }
  #### check how many zeros in Y
  if (sum(Y > 0) / length(Y) > 0.9) {
    warning("Too few zeros in the abundance data. The logistic component may be not accurate.")
  }
  if (sum(Y > 0) / length(Y) < 0.1) {
    warning("Too many zeros in the abundance data. The beta component may be not accurate.")
  }
  #### if the colnames are the same, jointly test the two component
  if (identical(colnames(logistic_cov), colnames(beta_cov))) {
    joint_test <- TRUE
  } else {
    joint_test <- FALSE
  }
  #### check if time_ind are the same for each subject_ind
  fit <- fit_zero_inflated_beta_random_effect(
    X = logistic_cov, Z = beta_cov, Y = Y,
    subject_ind = subject_ind, time_ind = time_ind, joint_test = joint_test, quad_n = quad_n
  )

  list(
    logistic_est_table = fit$logistic_est_table,
    logistic_s1_est = fit$logistic_s1_est,
    beta_est_table = fit$beta_est_table,
    beta_s2_est = fit$beta_s2_est,
    beta_v_est = fit$beta_v_est,
    loglikelihood = fit$loglikelihood,
    joint_p = fit$joint_p
  )
}
