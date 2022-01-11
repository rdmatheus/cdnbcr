test_that("rDNB under negative parameter information", {

  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, 0, 5)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, 0, 5), stats::runif(1, -5, 0)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, -5, 0)))), 100)
})


test_that("pgfDNB works", {

  s <- stats::runif(1, -1, 1)
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)

  expect_false(is.nan(pgfDNB(s, theta, phi, p, alpha)))

  # inconsistent arguments
  expect_true(is.nan(pgfDNB(s, -theta, phi, p, alpha)))
  expect_true(is.nan(pgfDNB(s, theta, -phi, p, alpha)))
  expect_true(is.nan(pgfDNB(s, theta, phi, -p, alpha)))
  expect_true(is.nan(pgfDNB(s, theta, phi, 2, alpha)))
  expect_true(is.nan(pgfDNB(s, theta, phi, p, -alpha)))
  expect_true(is.nan(pgfDNB(s, theta, phi, p, 2)))

  # Vectorization
  s <- matrix(stats::runif(6, -1, 1), ncol = 2)
  expect_equal(sum(pgfDNB(s, theta, phi, p, alpha) -  matrix(c(pgfDNB(s[1,1], theta, phi, p, alpha),
                                              pgfDNB(s[2,1], theta, phi, p, alpha),
                                              pgfDNB(s[3,1], theta, phi, p, alpha),
                                              pgfDNB(s[1,2], theta, phi, p, alpha),
                                              pgfDNB(s[2,2], theta, phi, p, alpha),
                                              pgfDNB(s[3,2], theta, phi, p, alpha)), ncol = 2)), 0)

  theta <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(pgfDNB(s, theta, phi, p, alpha) -  matrix(c(pgfDNB(s[1,1], theta[1,1], phi, p, alpha),
                                                  pgfDNB(s[2,1], theta[2,1], phi, p, alpha),
                                                  pgfDNB(s[3,1], theta[3,1], phi, p, alpha),
                                                  pgfDNB(s[1,2], theta[1,2], phi, p, alpha),
                                                  pgfDNB(s[2,2], theta[2,2], phi, p, alpha),
                                                  pgfDNB(s[3,2], theta[3,2], phi, p, alpha)), ncol = 2)), 0)


  p <- matrix(stats::runif(6, 0, 1), ncol = 2)
  expect_equal(sum(pgfDNB(s, theta, phi, p, alpha) -  matrix(c(pgfDNB(s[1,1], theta[1,1], phi, p[1,1], alpha),
                                                               pgfDNB(s[2,1], theta[2,1], phi, p[2,1], alpha),
                                                               pgfDNB(s[3,1], theta[3,1], phi, p[3,1], alpha),
                                                               pgfDNB(s[1,2], theta[1,2], phi, p[1,2], alpha),
                                                               pgfDNB(s[2,2], theta[2,2], phi, p[2,2], alpha),
                                                               pgfDNB(s[3,2], theta[3,2], phi, p[3,2], alpha)), ncol = 2)), 0)

  # Vectorization and inconsistent arguments
  theta[1, 2] <- -2
  theta[2, 1] <- -2
  p[1, 1] <- -1
  p[2, 2] <- 2
  expect_equal(sum(is.nan(pgfDNB(s, theta, phi, p, alpha))), 4)

})

test_that("pgf1DNB works", {

  s <- stats::runif(1, -1, 1)
  theta <- stats::runif(1, 0, 5)
  phi <- stats::runif(1, 0, 5)
  p <- stats::runif(1)
  alpha <- stats::runif(1)

  expect_false(is.nan(pgf1DNB(s, theta, phi, p, alpha)))

  # inconsistent arguments
  expect_true(is.nan(pgf1DNB(s, -theta, phi, p, alpha)))
  expect_true(is.nan(pgf1DNB(s, theta, -phi, p, alpha)))
  expect_true(is.nan(pgf1DNB(s, theta, phi, -p, alpha)))
  expect_true(is.nan(pgf1DNB(s, theta, phi, 2, alpha)))
  expect_true(is.nan(pgf1DNB(s, theta, phi, p, -alpha)))
  expect_true(is.nan(pgf1DNB(s, theta, phi, p, 2)))

  # Vectorization
  s <- matrix(stats::runif(6, -1, 1), ncol = 2)
  expect_equal(sum(pgf1DNB(s, theta, phi, p, alpha) -  matrix(c(pgf1DNB(s[1,1], theta, phi, p, alpha),
                                                               pgf1DNB(s[2,1], theta, phi, p, alpha),
                                                               pgf1DNB(s[3,1], theta, phi, p, alpha),
                                                               pgf1DNB(s[1,2], theta, phi, p, alpha),
                                                               pgf1DNB(s[2,2], theta, phi, p, alpha),
                                                               pgf1DNB(s[3,2], theta, phi, p, alpha)), ncol = 2)), 0)

  theta <- matrix(stats::runif(6, 0, 5), ncol = 2)
  expect_equal(sum(pgf1DNB(s, theta, phi, p, alpha) -  matrix(c(pgf1DNB(s[1,1], theta[1,1], phi, p, alpha),
                                                               pgf1DNB(s[2,1], theta[2,1], phi, p, alpha),
                                                               pgf1DNB(s[3,1], theta[3,1], phi, p, alpha),
                                                               pgf1DNB(s[1,2], theta[1,2], phi, p, alpha),
                                                               pgf1DNB(s[2,2], theta[2,2], phi, p, alpha),
                                                               pgf1DNB(s[3,2], theta[3,2], phi, p, alpha)), ncol = 2)), 0)


  p <- matrix(stats::runif(6, 0, 1), ncol = 2)
  expect_equal(sum(pgf1DNB(s, theta, phi, p, alpha) -  matrix(c(pgf1DNB(s[1,1], theta[1,1], phi, p[1,1], alpha),
                                                               pgf1DNB(s[2,1], theta[2,1], phi, p[2,1], alpha),
                                                               pgf1DNB(s[3,1], theta[3,1], phi, p[3,1], alpha),
                                                               pgf1DNB(s[1,2], theta[1,2], phi, p[1,2], alpha),
                                                               pgf1DNB(s[2,2], theta[2,2], phi, p[2,2], alpha),
                                                               pgf1DNB(s[3,2], theta[3,2], phi, p[3,2], alpha)), ncol = 2)), 0)

  # Vectorization and inconsistent arguments
  theta[1, 2] <- -2
  theta[2, 1] <- -2
  p[1, 1] <- -1
  p[2, 2] <- 2
  expect_equal(sum(is.nan(pgf1DNB(s, theta, phi, p, alpha))), 4)
})
