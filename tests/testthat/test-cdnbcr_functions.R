test_that("dcdnbcr works", {

  x <- 2
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)
  mu <- stats::runif(1, 0, 5)
  sigma <- stats::runif(1, 0, 5)

  expect_true(is.nan(dcdnbcr(-x, theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, -theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, -phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, -p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, -alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, alpha, -mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, alpha, mu, -sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, 2, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, 2, mu, sigma)))

  expect_error(dcdnbcr(x, c(2, 3, 4, 5), phi, p, alpha, mu, sigma))
  expect_error(dcdnbcr(x, theta, phi, c(0.2, 0.3, 0.4, 0.5), alpha, mu, sigma))
  expect_error(dcdnbcr(x, theta, phi, p, alpha, c(2, 3, 4, 5), sigma))

  x <- runif(10, -5, 5)
  expect_equal(sum(is.nan(dcdnbcr(x, theta, phi, p, alpha, mu, sigma))), sum(x < 0))

})


test_that("pcdnbcr works", {

  x <- 2
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)
  mu <- stats::runif(1, 0, 5)
  sigma <- stats::runif(1, 0, 5)

  expect_true(is.nan(pcdnbcr(-x, theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, -theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, -phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, -p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, -alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, alpha, -mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, alpha, mu, -sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, 2, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, 2, mu, sigma)))

  expect_error(pcdnbcr(x, c(2, 3, 4, 5), phi, p, alpha, mu, sigma))
  expect_error(pcdnbcr(x, theta, phi, c(0.2, 0.3, 0.4, 0.5), alpha, mu, sigma))
  expect_error(pcdnbcr(x, theta, phi, p, alpha, c(2, 3, 4, 5), sigma))

  x <- runif(10, -5, 5)
  expect_equal(sum(pcdnbcr(x, theta, phi, p, alpha, mu, sigma) == 0), sum(x < 0))

})


test_that("cure_rate works", {

  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)

  expect_false(is.nan(cure_rate(theta, phi, p, alpha)))
  expect_true(is.nan(cure_rate(-theta, phi, p, alpha)))
  expect_true(is.nan(cure_rate(theta, -phi, p, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, -p, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, 2, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, p, -alpha)))
  expect_true(is.nan(cure_rate(theta, phi, p, 2)))

  expect_equal(cure_rate(c(2, 3, 4), phi, p, alpha), c(cure_rate(2, phi, p, alpha),
                                                       cure_rate(3, phi, p, alpha),
                                                       cure_rate(4, phi, p, alpha)))
  expect_equal(cure_rate(theta, phi, c(0.2, 0.3, 0.4), alpha), c(cure_rate(theta, phi, 0.2, alpha),
                                                       cure_rate(theta, phi, 0.3, alpha),
                                                       cure_rate(theta, phi, 0.4, alpha)))

  theta <- runif(10, -5, 5)
  expect_equal(sum(is.nan(cure_rate(theta, phi, p, alpha))), sum(theta < 0))
  p <- runif(10, -5, 5)
  expect_equal(sum(is.nan(cure_rate(2, phi, p, alpha))), sum(p < 0 | p > 1))

})
