calc_beta_loglik <- function(para, Z.aug, Y, subject.n, time.n,
                            prod.mat = prod.mat,
                            Z.test.coeff.index,
                            gh.weights, gh.nodes, quad.n) {
  ### Y must be a vector here
  Y <- as.vector(Y)
  s2 <- para[1]
  v <- para[2]
  #### Z.test.coeff.index == TRUE means we want to test this parameter
  #### Therefore we want to keep it as 0 for LRT
  beta <- rep(NA, length(Z.test.coeff.index))
  beta[Z.test.coeff.index] <- 0 ## initialized as 0
  beta[!Z.test.coeff.index] <- para[-(1:2)][1:sum(!Z.test.coeff.index)]
  beta <- as.matrix(beta)
  u <- 1 / (1 + exp(-(Z.aug %*% beta[, rep(1, quad.n)] + gh.nodes * s2 * sqrt(2))))
  #### replace Y==0 with NA, so don't use them in the likelihood calculation
  # Y.tem[Y == 0] <- NA
  #### exp(log(A*B)) = exp(logA+logB)=A*B
  #### log likelihood
  #### first calculate the loglikelihood for the ith sbuject
  # log.i <- log(gamma(v)/(gamma(u*v)*gamma((1-u)*v))*Y^(u*v-1)*(1-Y)^((1-u)*v-1))
  suppressWarnings(
    log.i <- lgamma(v) - lgamma(u * v) - lgamma((1 - u) * v) + (u * v - 1) * log(Y) + ((1 - u) * v - 1) * log(1 - Y)
  )
  #######
  ## u -> 1 will cause warnings. lgamma((1-u)*v) = Inf.
  ## next line will hanle this Inf problem.
  #######
  # browser()
  #### logY==0 -> infinite, we want to exclude Y==0 in the likelihood calculation
  #### replace with infinite with 0, it will be ignored, e^0*e^log(p2)= p2
  log.i[is.infinite(log.i)] <- 0
  logL <- sum(
    log(rowSums(gh.weights / sqrt(pi) * exp(prod.mat %*% log.i),
      na.rm = TRUE
    )),
    na.rm = TRUE
  )

  -logL
}

#' Fit beta random effect
#'
#' @param Z FILL
#' @param Y FILL
#' @param subject.ind the subject index
#' @param time.ind the time index
#' @param quad.n number of points in gaussian quadrature
#' @param verbose a boolean to enable more output
#' @return a named list
#' \itemize{
#'   \item est.table
#'   \item s2.est
#'   \item v.est
#' }
#'
#' @importFrom stats nlminb pchisq
#' @importFrom statmod gauss.quad
fit_beta_random_effect <- function(Z = Z, Y = Y,
                                   subject.ind = subject.ind, time.ind = time.ind,
                                   quad.n = 30, verbose = FALSE) {
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  if (is.null(colnames(Z))) {
    colnames(Z) <- paste("var", seq_len(ncol(Z)), sep = "")
  }
  Z.aug <- cbind(intersept = 1, Z)

  est.table <- matrix(NA,
    ncol = 2, nrow = ncol(Z.aug),
    dimnames = list(colnames(Z.aug), c("Estimate", "Pvalue"))
  )

  subject.n <- length(unique(subject.ind))
  time.n <- length(unique(time.ind))
  prod.mat <- matrix(rep(c(rep(1, time.n), rep(0, subject.n * time.n)), subject.n)[1:(subject.n^2 * time.n)],
                     byrow = TRUE,
                     nrow = subject.n,
                     ncol = subject.n * time.n)

  #### generate quad points
  gherm <- gauss.quad(quad.n, kind = "hermite")
  gh.weights <- matrix(rep(gherm$weights, subject.n), nrow = subject.n, byrow = TRUE)
  gh.nodes <- matrix(rep(gherm$nodes, subject.n * time.n),
    nrow = subject.n * time.n, byrow = TRUE
  )

  #### re-order Z,Y so that values belongs to the same subject are together
  #### need the values in this format for loglikelihood calculation
  gind <- sort(subject.ind, index.return = TRUE)$ix
  Y <- Y[gind, ]
  Z <- Z[gind, ]
  Z.aug <- Z.aug[gind, ]
  ##### H1: estimate all parameters
  Z.test.coeff.index <- rep(FALSE, ncol(Z.aug))
  opt.H1 <- nlminb(
    start = c(1, 2, rep(0, sum(!Z.test.coeff.index))), ## s2,v,alpha
    objective = calc_beta_loglik,
    lower = c(
      0.00001, 0.00001,
      rep(-Inf, sum(!Z.test.coeff.index))
    ),
    upper = Inf,
    Z.test.coeff.index = Z.test.coeff.index,
    Y = Y, Z.aug = Z.aug, time.n = time.n, subject.n = subject.n,
    prod.mat = prod.mat,
    gh.weights = gh.weights, gh.nodes = gh.nodes, quad.n = quad.n,
    control = list(trace = ifelse(verbose, 2, 0))
  )
  #### save the estimatd results
  s2.est <- opt.H1$par[1]
  v.est <- opt.H1$par[2]
  beta.est <- opt.H1$par[-(1:2)]
  est.table[, "Estimate"] <- beta.est

  ####### H0
  for (test.i in seq_len(ncol(Z.aug))) {
    Z.test.coeff.index <- rep(FALSE, ncol(Z.aug))
    Z.test.coeff.index[test.i] <- TRUE
    opt.H0 <- nlminb(
      start = c(1, 2, rep(0, sum(!Z.test.coeff.index))), ## s2,alpha
      objective = calc_beta_loglik,
      lower = c(
        0.00001, 0.00001,
        rep(-Inf, sum(!Z.test.coeff.index))
      ),
      upper = c(Inf),
      Z.test.coeff.index = Z.test.coeff.index,
      Y = Y, Z.aug = Z.aug, time.n = time.n, subject.n = subject.n,
      prod.mat = prod.mat,
      gh.weights = gh.weights, gh.nodes = gh.nodes, quad.n = quad.n,
      control = list(trace = ifelse(verbose, 2, 0))
    )
    #### likelihood ratio test
    likelihodd.ratio <- -2 * (-opt.H0$objective - (-opt.H1$objective))
    LRT.p <- 1 - pchisq(likelihodd.ratio, df = 1)
    est.table[test.i, "Pvalue"] <- LRT.p
  }

  list(est.table = est.table, s2.est = s2.est, v.est = v.est)
}
