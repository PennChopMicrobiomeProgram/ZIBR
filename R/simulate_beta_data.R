#' Simulate beta data
#'
#' @param subject_n the number of subjects
#' @param time_n the number of time points
#' @param v FILL
#' @param beta FILL
#' @param Z FILL
#' @param s2 FILL
#' @param sim_seed the random seed with which to simulate the data
#' @return a named list
#' \itemize{
#'   \item Y
#'   \item Z
#'   \item c
#'   \item u
#'   \item v
#'   \item beta
#'   \item s2
#'   \item subject_ind
#'   \item time_ind
#' }
#'
#' @importFrom stats rnorm rbeta
simulate_beta_random_effect_data <- function(subject_n = 50, time_n = 5, v = 2,
                                             beta = as.matrix(c(-0.5, -0.5, 0.5)),
                                             Z = NA, s2 = 1, sim_seed = 100) {
  #### Beta regression
  if (length(NA) == 1 && any(is.na(Z))) {
    set.seed(sim_seed * 5000 + 10)
    Z <- as.matrix(data.frame(
      log.Time = as.matrix(log(rep(seq(1, time_n), subject_n))),
      Treatment = as.matrix(c(rep(0, subject_n * time_n / 2), rep(1, subject_n * time_n / 2)))
    ))
  }
  if (is.null(colnames(Z))) {
    Z <- as.matrix(Z)
    colnames(Z) <- paste("var", seq(1, ncol(Z)), sep = "")
  }

  set.seed(sim_seed * 5000 + 1)
  c <- as.matrix(rnorm(subject_n, mean = 0, sd = s2))
  c_rep <- as.matrix(as.vector(matrix(c, nrow = time_n, ncol = length(c), byrow = TRUE)))

  subject_ind <- as.vector(matrix(paste("Subject_", seq(1, subject_n), sep = ""),
                                  nrow = time_n,
                                  ncol = subject_n,
                                  byrow = TRUE))
  time_ind <- rep(seq(1, time_n), subject_n)

  z_aug <- cbind(intersept = 1, Z)

  logit_u <- z_aug %*% beta + c_rep

  u <- 1 / (1 + exp(-logit_u))
  ## v is the phi

  set.seed(sim_seed * 5000 + 4)
  Y <- rbeta(subject_n * time_n, shape1 = u * v, shape2 = (1 - u) * v)
  if (any(Y > 1 - 10^(-6))) {
    Y[Y > 1 - 10^(-6)] <- 1 - 10^(-6)
  }

  #### For test purpose
  #### betareg can not fit random effect model
  #### so set the s2 to a small value (small random effect)
  # library(betareg)
  # tdata <- data.frame(Y=Y,Z,SID=subject_ind)
  # gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='ML')
  # summary(gy)
  # gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BC')
  # summary(gy)
  # gy <- betareg(Y ~ log.Time + as.factor(Treatment) , data = tdata,type='BR')
  # summary(gy)
  #
  #
  # fit_beta_random_effect(Z=Z,Y=Y,subject_ind=subject_ind,time_ind=time_ind,
  #                      quad.n=30,verbose=FALSE)

  list(Y = Y, Z = Z, c = c, u = u, v = v, beta = beta, s2 = s2, subject_ind = subject_ind, time_ind = time_ind)
}
