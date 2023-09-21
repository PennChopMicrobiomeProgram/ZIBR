#' @importFrom stats rnorm rbinom
simulate_logistic_data <- function() {
  #### logistic regression
  sim.seed <- 18642
  time.n <- 5
  subject.n <- 20
  s1 <- 0.5
  alpha <- as.matrix(c(0, 0.5, -1))
  ######
  set.seed(sim.seed + 10)
  X <- as.matrix(data.frame(
    log.Time = as.matrix(log(rep(seq(1, time.n), subject.n))),
    Treatment = as.matrix(c(rep(0, subject.n * time.n / 2), rep(1, subject.n * time.n / 2)))
  ))
  ######
  set.seed(sim.seed + 2)
  b <- as.matrix(rnorm(subject.n, mean = 0, sd = s1))
  b.rep <- as.matrix(as.vector(matrix(b, nrow = time.n, ncol = length(b), byrow = TRUE)))
  #####
  subject.ind <- as.vector(matrix(paste("Subject_", seq(1, subject.n), sep = ""), nrow = time.n, ncol = subject.n, byrow = TRUE))
  time.ind <- rep(seq(1, time.n), subject.n)
  ######
  X.aug <- cbind(interespt = 1, X)
  ######
  logit.p <- X.aug %*% alpha + b.rep
  ######
  p <- 1 / (1 + exp(-logit.p))
  ######
  set.seed(sim.seed + 3)
  Y <- rbinom(subject.n * time.n, 1, p)

  list(X = X, Y = Y, b = b, subject.ind = subject.ind, time.ind = time.ind)
}
