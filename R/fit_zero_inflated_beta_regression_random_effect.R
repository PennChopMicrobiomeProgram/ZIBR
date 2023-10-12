calc_zibeta_loglik <- function(para,
                              subject_n, time_n,
                              X.aug, Z.aug, Y,
                              X.test.coeff.index,
                              Z.test.coeff.index,
                              prod.mat,
                              gh.weights, gh.nodes,
                              quad.n = quad.n) {
  #### s1, s2, v will always be estimated
  logistic.para <- para[1:(sum(!X.test.coeff.index) + 1)] ## s1, alpha
  beta.para <- para[-(1:(sum(!X.test.coeff.index) + 1))] ## s2,v, beta
  logistic.logLike <- calc_logistic_loglik(
    para = logistic.para,
    X.aug = X.aug, Y = Y,
    subject.n = subject_n,
    time.n = time_n,
    prod.mat = prod.mat,
    X.test.coeff.index = X.test.coeff.index,
    gh.weights = gh.weights, gh.nodes = gh.nodes,
    quad.n = quad.n
  )
  beta.logLike <- calc_beta_loglik(
    para = beta.para,
    Z.aug = Z.aug, Y = Y,
    subject.n = subject_n,
    time.n = time_n,
    Z.test.coeff.index = Z.test.coeff.index,
    prod.mat = prod.mat,
    gh.weights = gh.weights, gh.nodes = gh.nodes,
    quad.n = quad.n
  )

  logistic.logLike + beta.logLike
}

#' Fit zero inflated beta random effect
#'
#' @param X FILL
#' @param Z FILL
#' @param Y FILL
#' @param subject_ind the subject index
#' @param time_ind the time index
#' @param component_wise_test boolean to run component-wise test
#' @param joint_test boolean to run joint test
#' @param quad_n number of points in gaussian quadrature
#' @param verbose a boolean to enable more output
#' @return a named list
#' \itemize{
#'   \item logistic_est_table
#'   \item logistic_s1_est
#'   \item beta_est_table
#'   \item beta_s2_est
#'   \item beta_v_est
#'   \item loglikelihood
#'   \item joint_p
#' }
#'
#' @importFrom stats nlminb pchisq
#' @importFrom statmod gauss.quad
fit_zero_inflated_beta_random_effect <- function(X = X, Z = Z, Y = Y,
                                                 subject_ind = subject_ind, time_ind = time_ind,
                                                 component_wise_test = TRUE,
                                                 joint_test = TRUE,
                                                 quad_n = 30, verbose = FALSE) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  if (is.null(colnames(X))) {
    colnames(X) <- paste("var", seq_len(ncol(X)), sep = "")
  }
  if (is.null(colnames(Z))) {
    colnames(Z) <- paste("var", seq_len(ncol(Z)), sep = "")
  }
  X.aug <- cbind(intersept = 1, X)
  Z.aug <- cbind(intersept = 1, Z)

  subject.n <- length(unique(subject_ind))
  time.n <- length(unique(time_ind))
  prod.mat <- matrix(rep(c(rep(1, time.n), rep(0, subject.n * time.n)), subject.n)[1:(subject.n^2 * time.n)],
                     byrow = TRUE,
                     nrow = subject.n,
                     ncol = subject.n * time.n)

  #### generate quad points
  gherm <- gauss.quad(quad_n, kind = "hermite")
  gh.weights <- matrix(rep(gherm$weights, subject.n), nrow = subject.n, byrow = TRUE)
  gh.nodes <- matrix(rep(gherm$nodes, subject.n * time.n),
    nrow = subject.n * time.n, byrow = TRUE
  )

  ###### estimate and test each parameter
  if (component_wise_test) {
    logistic.fit <- fit_logistic_random_effect(
      X = X, Y = Y,
      subject.ind = subject_ind, time.ind = time_ind,
      quad.n = quad_n, verbose = verbose
    )
    beta.fit <- fit_beta_random_effect(
      Z = Z, Y = Y, subject.ind = subject_ind, time.ind = time_ind,
      quad.n = quad_n, verbose = verbose
    )
  }

  ##### jointly test each variable in logistic and beta component
  if (joint_test) {
    ##### H1
    X.test.coeff.index <- rep(FALSE, ncol(X.aug))
    Z.test.coeff.index <- rep(FALSE, ncol(Z.aug))
    opt.H1 <- nlminb(
      start = c(
        c(1, rep(0, sum(!X.test.coeff.index))), ## s1,alpha
        c(1, 2, rep(0, sum(!Z.test.coeff.index)))
      ), ## s2,v,beta
      objective = calc_zibeta_loglik,
      lower = c(
        c(0.00001, rep(-Inf, sum(!X.test.coeff.index))),
        c(0.00001, 0.00001, rep(-Inf, sum(!Z.test.coeff.index)))
      ),
      upper = Inf,
      X.test.coeff.index = X.test.coeff.index,
      Z.test.coeff.index = Z.test.coeff.index,
      Y = Y, X.aug = X.aug, Z.aug = Z.aug, time_n = time.n, subject_n = subject.n,
      prod.mat = prod.mat,
      gh.weights = gh.weights, gh.nodes = gh.nodes,
      quad.n = quad_n,
      control = list(trace = ifelse(verbose, 2, 0))
    )
    ###### H0:set corresponding regression coefficients to zero, excluding intercept
    joint_p <- rep(NA, ncol(X.aug))
    for (i in 2:ncol(X.aug)) {
      X.test.coeff.index <- rep(FALSE, ncol(X.aug))
      Z.test.coeff.index <- rep(FALSE, ncol(Z.aug))
      X.test.coeff.index[i] <- TRUE
      Z.test.coeff.index[i] <- TRUE
      opt.H0 <- nlminb(
        start = c(
          c(1, rep(0, sum(!X.test.coeff.index))), ## s1,alpha
          c(1, 2, rep(0, sum(!Z.test.coeff.index)))
        ), ## s2,v,beta
        objective = calc_zibeta_loglik,
        lower = c(
          c(0.00001, rep(-Inf, sum(!X.test.coeff.index))),
          c(0.00001, 0.00001, rep(-Inf, sum(!Z.test.coeff.index)))
        ),
        upper = Inf,
        X.test.coeff.index = X.test.coeff.index,
        Z.test.coeff.index = Z.test.coeff.index,
        Y = Y, X.aug = X.aug, Z.aug = Z.aug, time_n = time.n, subject_n = subject.n,
        prod.mat = prod.mat,
        gh.weights = gh.weights, gh.nodes = gh.nodes,
        quad.n = quad_n,
        control = list(trace = ifelse(verbose, 2, 0))
      )
      likelihodd.ratio <- -2 * (-opt.H0$objective - (-opt.H1$objective))
      joint_p[i] <- 1 - pchisq(likelihodd.ratio, df = sum(X.test.coeff.index) + sum(Z.test.coeff.index))
    }
    names(joint_p) <- c("overall", colnames(X))
    ## remove the 'overall'
    joint_p <- joint_p[-1]
  }

  return_values <- list(logistic_est_table = NA,
                        logistic_s1_est = NA,
                        beta_est_table = NA,
                        beta_s2_est = NA,
                        beta_v_est = NA,
                        loglikelihood = NA,
                        joint_p = NA)

  if (component_wise_test) {
    return_values[["logistic_est_table"]] <- logistic.fit$est.table
    return_values[["logistic_s1_est"]] <- logistic.fit$s1.est
    return_values[["beta_est_table"]] <- beta.fit$est.table
    return_values[["beta_s2_est"]] <- beta.fit$s2.est
    return_values[["beta_v_est"]] <- beta.fit$v.est
  }
  if (joint_test) {
    return_values[["loglikelihood"]] <- -opt.H1$objective
    return_values[["joint_p"]] <- joint_p
  }

  return_values
}
