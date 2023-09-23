#' Simulate logistic data
#'
#' @return a named list
#' \itemize{
#'   \item X
#'   \item Y
#'   \item b
#'   \item subject_ind
#'   \item time_ind
#' }
#'
#' @importFrom stats rnorm rbinom
simulate_logistic_data <- function() {
  #### logistic regression
  sim_seed <- 18642
  time_n <- 5
  subject_n <- 20
  s1 <- 0.5
  alpha <- as.matrix(c(0, 0.5, -1))

  set.seed(sim_seed + 10)
  X <- as.matrix(data.frame(
    log.Time = as.matrix(log(rep(seq(1, time_n), subject_n))),
    Treatment = as.matrix(c(rep(0, subject_n * time_n / 2), rep(1, subject_n * time_n / 2)))
  ))

  set.seed(sim_seed + 2)
  b <- as.matrix(rnorm(subject_n, mean = 0, sd = s1))
  b_rep <- as.matrix(as.vector(matrix(b, nrow = time_n, ncol = length(b), byrow = TRUE)))

  subject_ind <- as.vector(matrix(paste("Subject_", seq(1, subject_n), sep = ""),
                                  nrow = time_n,
                                  ncol = subject_n,
                                  byrow = TRUE))
  time_ind <- rep(seq(1, time_n), subject_n)

  x_aug <- cbind(interespt = 1, X)

  logit_p <- x_aug %*% alpha + b_rep

  p <- 1 / (1 + exp(-logit_p))

  set.seed(sim_seed + 3)
  Y <- rbinom(subject_n * time_n, 1, p)

  list(X = X, Y = Y, b = b, subject_ind = subject_ind, time_ind = time_ind)
}
