real_expected <- list(logistic.est.table = structure(c(2.69732340911363, 0.200018071642455, 6.17105397782147e-06, 0.887710710500198), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "Treatment"), c("Estimate", "Pvalue"))), logistic.s1.est = 3.27876357265477, beta.est.table = structure(c(-2.78605402671509, -0.32334269957996, 0, 0.206300170582503), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "Treatment"), c("Estimate", "Pvalue"))), beta.s2.est = 0.538923350735293, beta.v.est = 7.62238420518043, loglikelihood = 321.901751058239, joint.p = c(Treatment = 0.445494681803514))

test_that("zibr main function works on package data", {
  zibr.fit <- zibr(logistic.cov = data.frame(Treatment=ibd$Treatment),
                   beta.cov = data.frame(Treatment=ibd$Treatment),
                   Y = ibd$Abundance, subject.ind = ibd$Subject,
                   time.ind = ibd$Time)

  expect_equal(names(zibr.fit), c("logistic.est.table", "logistic.s1.est", "beta.est.table", "beta.s2.est", "beta.v.est", "loglikelihood", "joint.p"))
  expect_equal(zibr.fit, real_expected, tolerance = 1e-3)
})

set.seed(19683)

sim <- simulate_zero_inflated_beta_random_effect_data(
  subject.n=100,time.n=5,
  X = as.matrix(c(rep(0,50*5),rep(1,50*5))),
  Z = as.matrix(c(rep(0,50*5),rep(1,50*5))),
  alpha = as.matrix(c(-0.5,1)),
  beta = as.matrix(c(-0.5,0.5)),
  s1 = 1,s2 = 0.8,
  v = 5,
  sim.seed=100)

sim_expected <- list(logistic.est.table = structure(c(-0.571900346260645, 0.828072571379781, 0.0069577928092327, 0.00547006539016526), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "var1"), c("Estimate", "Pvalue"))), logistic.s1.est = 1.06801391077711, beta.est.table = structure(c(-0.593090577296489, 0.591745582740125, 4.18162741963046e-05, 0.00169945318503251), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "var1"), c("Estimate", "Pvalue"))), beta.s2.est = 0.630386186272789, beta.v.est = 4.73406414312415, loglikelihood = NA, joint.p = NA)

sim_expected_cov <- list(logistic.est.table = structure(c(-0.571900346260645, 0.828072571379781, 0.0069577928092327, 0.00547006539016526), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "var1"), c("Estimate", "Pvalue"))), logistic.s1.est = 1.06801391077711, beta.est.table = structure(c(-0.593090577296489, 0.591745582740125, 4.18162741963046e-05, 0.00169945318503251), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "var1"), c("Estimate", "Pvalue"))), beta.s2.est = 0.630386186272789, beta.v.est = 4.73406414312415, loglikelihood = -286.585518147651, joint.p = c(var1 = 0.000153328547781939))

test_that("zibr main function works on simulated data", {
  zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$Z, Y = sim$Y, subject.ind = sim$subject.ind,time.ind = sim$time.ind)

  expect_equal(zibr.fit, sim_expected, tolerance = 1e-3)
})

test_that("zibr main function works on simulated data with same covariates", {
  zibr.fit <- zibr(logistic.cov = sim$X, beta.cov = sim$X, Y = sim$Y, subject.ind = sim$subject.ind,time.ind = sim$time.ind)

  expect_equal(zibr.fit, sim_expected_cov, tolerance = 1e-3)
})
