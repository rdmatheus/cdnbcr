test_that("dcbern works", {

  expect_equal(dcbern(stats::runif(1, -10, 10), stats::runif(1), stats::runif(1)), 0)
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), -stats::runif(1), stats::runif(1))))
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), stats::runif(1), -stats::runif(1))))
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), stats::runif(1, -1, 0), stats::runif(1))))
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), stats::runif(1, 1, 2), stats::runif(1))))
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), stats::runif(1), stats::runif(1, -1, 0))))
  expect_true(is.nan(dcbern(stats::runif(1, -10, 10), stats::runif(1), stats::runif(1, 1, 2))))
})
