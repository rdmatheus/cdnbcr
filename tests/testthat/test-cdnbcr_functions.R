test_that("dcdnbcr works", {

  x <- 2
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)
  mu <- stats::runif(1, 0, 5)
  sigma <- stats::runif(1, 0, 5)

  # Inconsistent arguments
  expect_true(is.nan(dcdnbcr(-x, theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, -theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, -phi, p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, -p, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, -alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, alpha, -mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, alpha, mu, -sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, 2, alpha, mu, sigma)))
  expect_true(is.nan(dcdnbcr(x, theta, phi, p, 2, mu, sigma)))

  # Vectorization
  x <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(dcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(dcdnbcr(x[1,1], theta, phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[2,1], theta, phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[3,1], theta, phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[1,2], theta, phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[2,2], theta, phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[3,2], theta, phi, p, alpha, mu, sigma)), ncol = 2)), 0)

  theta <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(dcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(dcdnbcr(x[1,1], theta[1,1], phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[2,1], theta[2,1], phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[3,1], theta[3,1], phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[1,2], theta[1,2], phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[2,2], theta[2,2], phi, p, alpha, mu, sigma),
                                                                           dcdnbcr(x[3,2], theta[3,2], phi, p, alpha, mu, sigma)), ncol = 2)), 0)


  p <- matrix(stats::runif(6, 0, 1), ncol = 2)
  expect_equal(sum(dcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(dcdnbcr(x[1,1], theta[1,1], phi, p[1,1], alpha, mu, sigma),
                                                                           dcdnbcr(x[2,1], theta[2,1], phi, p[2,1], alpha, mu, sigma),
                                                                           dcdnbcr(x[3,1], theta[3,1], phi, p[3,1], alpha, mu, sigma),
                                                                           dcdnbcr(x[1,2], theta[1,2], phi, p[1,2], alpha, mu, sigma),
                                                                           dcdnbcr(x[2,2], theta[2,2], phi, p[2,2], alpha, mu, sigma),
                                                                           dcdnbcr(x[3,2], theta[3,2], phi, p[3,2], alpha, mu, sigma)), ncol = 2)), 0)

  theta[1, 2] <- -2
  theta[2, 1] <- -2
  p[1, 1] <- -1
  p[2, 2] <- 2
  expect_equal(sum(is.nan(dcdnbcr(x, theta, phi, p, alpha, mu, sigma))), 4)

  x <- runif(10, -5, 5)
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)
  mu <- stats::runif(1, 0, 5)
  sigma <- stats::runif(1, 0, 5)
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

  # Inconsistent arguments
  expect_true(is.nan(pcdnbcr(-x, theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, -theta, phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, -phi, p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, -p, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, -alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, alpha, -mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, alpha, mu, -sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, 2, alpha, mu, sigma)))
  expect_true(is.nan(pcdnbcr(x, theta, phi, p, 2, mu, sigma)))

  # Vectorization
  x <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(pcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(pcdnbcr(x[1,1], theta, phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[2,1], theta, phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[3,1], theta, phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[1,2], theta, phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[2,2], theta, phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[3,2], theta, phi, p, alpha, mu, sigma)), ncol = 2)), 0)

  theta <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(pcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(pcdnbcr(x[1,1], theta[1,1], phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[2,1], theta[2,1], phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[3,1], theta[3,1], phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[1,2], theta[1,2], phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[2,2], theta[2,2], phi, p, alpha, mu, sigma),
                                                                           pcdnbcr(x[3,2], theta[3,2], phi, p, alpha, mu, sigma)), ncol = 2)), 0)


  p <- matrix(stats::runif(6, 0, 1), ncol = 2)
  expect_equal(sum(pcdnbcr(x, theta, phi, p, alpha, mu, sigma) -  matrix(c(pcdnbcr(x[1,1], theta[1,1], phi, p[1,1], alpha, mu, sigma),
                                                                           pcdnbcr(x[2,1], theta[2,1], phi, p[2,1], alpha, mu, sigma),
                                                                           pcdnbcr(x[3,1], theta[3,1], phi, p[3,1], alpha, mu, sigma),
                                                                           pcdnbcr(x[1,2], theta[1,2], phi, p[1,2], alpha, mu, sigma),
                                                                           pcdnbcr(x[2,2], theta[2,2], phi, p[2,2], alpha, mu, sigma),
                                                                           pcdnbcr(x[3,2], theta[3,2], phi, p[3,2], alpha, mu, sigma)), ncol = 2)), 0)


  theta[1, 2] <- -2
  theta[2, 1] <- -2
  p[1, 1] <- -1
  p[2, 2] <- 2
  expect_equal(sum(is.nan(pcdnbcr(x, theta, phi, p, alpha, mu, sigma))), 4)

  x <- runif(10, -5, 5)
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)
  mu <- stats::runif(1, 0, 5)
  sigma <- stats::runif(1, 0, 5)
  expect_equal(sum(pcdnbcr(x, theta, phi, p, alpha, mu, sigma) == 0), sum(x < 0))

})


test_that("cure_rate works", {

  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)

  expect_false(is.nan(cure_rate(theta, phi, p, alpha)))

  # Inconsistent arguments
  expect_true(is.nan(cure_rate(-theta, phi, p, alpha)))
  expect_true(is.nan(cure_rate(theta, -phi, p, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, -p, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, 2, alpha)))
  expect_true(is.nan(cure_rate(theta, phi, p, -alpha)))
  expect_true(is.nan(cure_rate(theta, phi, p, 2)))

  # Vectorization
  expect_equal(cure_rate(c(2, 3, 4), phi, p, alpha), c(cure_rate(2, phi, p, alpha),
                                                       cure_rate(3, phi, p, alpha),
                                                       cure_rate(4, phi, p, alpha)))
  expect_equal(cure_rate(theta, phi, c(0.2, 0.3, 0.4), alpha), c(cure_rate(theta, phi, 0.2, alpha),
                                                       cure_rate(theta, phi, 0.3, alpha),
                                                       cure_rate(theta, phi, 0.4, alpha)))

  # Inconsistent arguments and vectorization
  theta <- runif(10, -5, 5)
  expect_equal(sum(is.nan(cure_rate(theta, phi, p, alpha))), sum(theta < 0))
  p <- runif(10, -5, 5)
  expect_equal(sum(is.nan(cure_rate(2, phi, p, alpha))), sum(p < 0 | p > 1))

})
