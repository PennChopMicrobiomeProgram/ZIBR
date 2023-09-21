#' Simulate data according to zero-inflated beta random effects model
#'
#' @param subject_n number of subjects
#' @param time_n number of time points for each subject
#' @param v the dispersion parameter in beta component
#' @param alpha the coefficients in logistic component
#' @param beta the coefficients in beta component
#' @param X the covariates in logistic component
#' @param Z the covariates in beta component
#' @param s1 the stardard deviation of random effect in logistic component
#' @param s2 the stardard deviation of random effect in beta component
#' @param sim_seed the random seed
#' @return Y the bacterial abundance generated from the model
#' @return X the covariates in logistic component
#' @return Z the covariates in beta component
#' @return alpha the coefficients in logistic component
#' @return beta the coefficients in beta component
#' @return s1 the stardard deviation of random effect in logistic component
#' @return s2 the stardard deviation of random effect in beta component
#' @return subject_ind the IDs for each subject
#' @return time_ind time points
#' @importFrom stats rnorm rbinom rbeta
#' @export
#' @examples
#' \dontrun{
#' simulate_zero_inflated_beta_random_effect_data(
#'   subject_n = 100, time_n = 5,
#'   X = as.matrix(c(rep(0, 50 * 5), rep(1, 50 * 5))),
#'   alpha = as.matrix(c(-0.5, 1)),
#'   beta = as.matrix(c(-0.5, 0.5)),
#'   s1 = 1, s2 = 0.8,
#'   v = 5,
#'   sim_seed = 100
#' )
#' }
simulate_zero_inflated_beta_random_effect_data <- function(subject_n = 50, time_n = 5, v = 2,
                                                           alpha = as.matrix(c(0, 0.5, -1)),
                                                           beta = as.matrix(c(-0.5, -0.5, 0.5)),
                                                           X = NA, Z = NA,
                                                           s1 = 0.2, s2 = 0.2, sim_seed = 100) {
  ######
  if (length(NA) == 1 & any(is.na(X))) {
    set.seed(sim_seed * 5000 + 10)
    X <- as.matrix(data.frame(
      log.Time = as.matrix(log(rep(seq(1, time_n), subject_n))),
      Treatment = as.matrix(c(rep(0, subject_n * time_n / 2), rep(1, subject_n * time_n / 2)))
    ))
    # X <- as.matrix(runif(subject_n*time_n,-1,1))
    # X  <- as.matrix(c(rep(0,N*time_n/2),rep(1,N*time_n/2)))
  }
  if (is.null(colnames(X))) {
    X <- as.matrix(X)
    colnames(X) <- paste("var", seq(1, ncol(X)), sep = "")
  }
  if (length(NA) == 1 & any(is.na(Z))) {
    Z <- X
  }


  set.seed(sim_seed * 5000 + 1)
  b <- as.matrix(rnorm(subject_n, mean = 0, sd = s1))
  b_rep <- as.matrix(as.vector(matrix(b, nrow = time_n, ncol = length(b), byrow = TRUE)))
  set.seed(sim_seed * 5000 + 2)
  c <- as.matrix(rnorm(subject_n, mean = 0, sd = s2))
  c_rep <- as.matrix(as.vector(matrix(c, nrow = time_n, ncol = length(c), byrow = TRUE)))
  #####
  subject_ind <- as.vector(matrix(paste("Subject_", seq(1, subject_n), sep = ""), nrow = time_n, ncol = subject_n, byrow = TRUE))
  time_ind <- rep(seq(1, time_n), subject_n)
  ######
  x_aug <- cbind(interespt = 1, X)
  z_aug <- cbind(intersept = 1, Z)
  # browser()
  logit_p <- X_aug %*% alpha + b_rep
  logit_u <- z_aug %*% beta + c_rep
  # print(z_aug %*% beta)
  ##### beta can not be too big, otherwise u will be close to 1

  ######
  p <- 1 / (1 + exp(-logit_p))
  u <- 1 / (1 + exp(-logit_u))
  # set.seed(sim_seed+2)
  # v <- runif(1,0,5) ## v is the phi
  if (is.na(v)) {
    v <- 2
  }


  ######
  Y <- rep(NA, subject_n * time_n)
  set.seed(sim_seed * 5000 + 3)
  ind_mix <- rbinom(length(Y), 1, p)
  for (i in 1:length(Y)) {
    if (ind_mix[i] == 0) {
      Y[i] <- 0
    }
    if (ind_mix[i] == 1) {
      # rbeta(n, shape1, shape2, ncp = 0)
      set.seed(sim_seed * 5000 + i)
      Y[i] <- rbeta(1, shape1 = u[i] * v, shape2 = (1 - u[i]) * v)
      if (Y[i] > 1 - 10^(-6)) {
        Y[i] <- 1 - 10^(-6)
      }
      # while(round(Y[i],6)==1){
      #  Y[i] = rbeta(1, shape1 = u[i]*v, shape2=(1-u[i])*v)
      # }
    }
  }

  list(Y = Y, X = X, Z = Z, b = b, c = c, u = u, v = v, alpha = alpha, beta = beta, s1 = s1, s2 = s2, subject_ind = subject_ind, time_ind = time_ind)
}
