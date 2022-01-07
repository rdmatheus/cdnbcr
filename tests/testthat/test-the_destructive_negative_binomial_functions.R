test_that("rDNB works", {

  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, 0, 5)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, 0, 5), stats::runif(1, -5, 0)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, -5, 0)))), 100)
})
