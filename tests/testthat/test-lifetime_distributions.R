test_that("correct specification of the reparameterized Weibull distribution", {

  mu <- runif(1, 0, 10)
  sigma <- runif(1, 0, 10)
  y <- rrweibull(10000, mu, sigma)
  expect_lt(median(ecdf(y)(y) - prweibull(y, mu, sigma)), 0.05)

  # par(mfrow = c(1, 2))
  # hist(y, prob = TRUE)
  # curve(drweibull(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
  # plot(ecdf(y))
  # curve(prweibull(x, mu, sigma), add = TRUE, col = 2, lwd = 2)
  # par(mfrow = c(1, 1))

})
