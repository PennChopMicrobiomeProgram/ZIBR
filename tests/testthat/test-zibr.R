E <- list(logistic.est.table = structure(c(2.69732340911363, 0.200018071642455, 6.17105397782147e-06, 0.887710710500198), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "Treatment"), c("Estimate", "Pvalue"))), logistic.s1.est = 3.27876357265477, beta.est.table = structure(c(-2.78605402671509, -0.32334269957996, 0, 0.206300170582503), .Dim = c(2L, 2L), .Dimnames = list(c("intersept", "Treatment"), c("Estimate", "Pvalue"))), beta.s2.est = 0.538923350735293, beta.v.est = 7.62238420518043, loglikelihood = 321.901751058239, joint.p = c(Treatment = 0.445494681803514))

test_that("zibr main function works", {
  zibr.fit <- zibr(logistic.cov = data.frame(Treatment=ibd.data$Treatment),
                   beta.cov = data.frame(Treatment=ibd.data$Treatment),
                   Y = ibd.data$Abundance, subject.ind = ibd.data$Subject,
                   time.ind = ibd.data$Time)
  print("HEYS")
  dput(zibr.fit)

  expect_equal(names(zibr.fit), c("logistic.est.table", "logistic.s1.est", "beta.est.table", "beta.s2.est", "beta.v.est", "loglikelihood", "joint.p"))
  expect_equal(zibr.fit, E)
})