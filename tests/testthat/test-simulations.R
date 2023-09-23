expected_bre_fit <- list(
  est.table = structure(
    c(-0.60955243747617, -0.443501946207444, 0.736943004541776,
      0.0104576873063189, 0.000414042620671129, 0.0123799088684966),
    .Dim = 3:2,
    .Dimnames = list(
      c("intersept", "log.Time", "Treatment"),
      c("Estimate", "Pvalue")
    )
  ),
  s2.est = 0.881964535293812,
  v.est = 2.14015032324312
)

test_that("simulate_beta_data works", {
  beta_data <- simulate_beta_random_effect_data()

  # betareg can not fit random effect model so set the s2 to a small value (small random effect)
  tdata <- data.frame(Y = beta_data$Y, beta_data$Z, SID = beta_data$subject_ind)
  gy <- betareg::betareg(beta_data$Y ~ log.Time + as.factor(Treatment), data = tdata, type = "ML")
  #print(summary(gy))
  gy <- betareg::betareg(beta_data$Y ~ log.Time + as.factor(Treatment), data = tdata, type = "BC")
  #print(summary(gy))
  gy <- betareg::betareg(beta_data$Y ~ log.Time + as.factor(Treatment), data = tdata, type = "BR")
  #print(summary(gy))

  res <- fit_beta_random_effect(Z = beta_data$Z,
                                Y = beta_data$Y,
                                subject.ind = beta_data$subject_ind,
                                time.ind = beta_data$time_ind,
                                quad.n = 30,
                                verbose = FALSE)

  expect_equal(res, expected_bre_fit, tolerance = 1e-3)
})

expected_lre_fit <- list(
  est.table = structure(
    c(0.306795205483351, 0.280454254105502,
      -0.667159092633291, 0.431194361789462, 0.266160176079596, 0.123214729652536),
    .Dim = 3:2,
    .Dimnames = list(c("intersept", "log.Time", "Treatment"), c("Estimate", "Pvalue"))
  ),
  s1.est = 1.12114408844058
)

test_that("simulate_logistic_data works", {
  logistic_data <- simulate_logistic_data()

  tdata <- data.frame(Y = logistic_data$Y, logistic_data$X, SID = logistic_data$subject_ind)
  lme.fit <- lme4::glmer(as.factor(logistic_data$Y) ~ log.Time + as.factor(Treatment) + (1 | SID),
                         data = tdata,
                         family = binomial,
                         control = lme4::glmerControl(optimizer = "bobyqa"),
                         nAGQ = 10)
  summary(lme.fit)

  res <- fit_logistic_random_effect(X = logistic_data$X,
                                    Y = logistic_data$Y,
                                    subject.ind = logistic_data$subject_ind,
                                    time.ind = logistic_data$time_ind)
  dput(res)

  expect_equal(res, expected_lre_fit, tolerance = 1e-3)
})
