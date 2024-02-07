test_that(desc = "basic chronospace output", code = {
  data("echinoid_dates")
  cspace <- chronospace(data_ages = echinoid_dates)

  expect_match(class(cspace), "chronospace")
  expect_equal(length(cspace), ncol(echinoid_dates$factors) + 1)
  expect_true("Total_vartable" %in% names(cspace))
})


test_that(desc = "stored factor information: ordination", code = {
  data("echinoid_dates")
  cspace <- chronospace(data_ages = echinoid_dates)

  #dimensions
  expect_true(nrow(cspace$factor_A$ordination$x) == nrow(echinoid_dates$ages))
  expect_true(ncol(cspace$factor_A$ordination$x) == nlevels(echinoid_dates$factors$factor_A) - 1)

  #scores, eigenvecs and eigenvals
  gx <- apply(cspace$factor_A$data$ages, 2, FUN = tapply, INDEX = cspace$factor_A$data$groups, mean)
  bgpca <- stats::prcomp(gx)

  x <- cspace$factor_A$ordination$x
  rotation <- cspace$factor_A$ordination$rotation
  values <- cspace$factor_A$ordination$values

  expect_true(all(abs(round(apply(x, 2, FUN = tapply, INDEX = cspace$factor_A$data$groups, mean)[,1:2],5)) == abs(round(bgpca$x[,1:2],5))))
  expect_true(all(round(abs(bgpca$rotation[,1:2]),5) == round(abs(rotation[,1:2]),5)))
  expect_true(all(round(values[1:2],5) == round(bgpca$sdev[1:2]^2, 5)))

})


test_that(desc = "stored factor information: sum of squares", code = {
  data("echinoid_dates")
  cspace <- chronospace(data_ages = echinoid_dates)

  #sum of squares (first bgPC only)
  Y <- cspace$factor_A$ordination$x[,1]
  A <- echinoid_dates$factors$factor_A
  B <- echinoid_dates$factors$factor_B

  nadj <- length(Y) - 1
  totvar <- cspace$factor_A$ssq$totvar

  X <- model.matrix(~ A + B)
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  fitted <- X %*% beta
  residuals <- Y - fitted
  SSr_perc <- 100 * ((t(residuals) %*% residuals) / nadj) / totvar

  SSa_b <- t(fitted) %*% fitted
  X <- model.matrix(~ B)
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  fitted <- X %*% beta
  SSb0 <- t(fitted) %*% fitted
  SSa_perc <- 100 * ((SSa_b - SSb0) / nadj) / totvar

  X <- model.matrix(~ A)
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  fitted <- X %*% beta
  SSa0 <- t(fitted) %*% fitted
  SSb_perc <- 100 * ((SSa_b - SSa0) / nadj) / totvar

  expect_true(all(round(c(SSa_perc, SSb_perc, SSr_perc),5) == round(as.numeric(cspace$factor_A$ssq$vartable[1,2:4]),5)))

})



