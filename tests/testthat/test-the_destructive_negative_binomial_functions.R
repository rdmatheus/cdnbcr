test_that("rDNB under negative parameter information", {

  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, 0, 5)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, 0, 5), stats::runif(1, -5, 0)))), 100)
  expect_equal(sum(is.nan(rDNB(100, stats::runif(1, -5, 0), stats::runif(1, -5, 0)))), 100)
})


test_that("pgfDNB works", {

  expect_false(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), -stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), -stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1, -1, 0), stats::runif(1))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1, 1, 2), stats::runif(1))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1, -1, 0))))
  expect_true(is.nan(pgfDNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1, 1, 2))))
})

test_that("pgf1DNB works", {

  expect_false(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), -stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), -stats::runif(1, 0, 5), stats::runif(1), stats::runif(1))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1, -1, 0), stats::runif(1))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1, 1, 2), stats::runif(1))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1, -1, 0))))
  expect_true(is.nan(pgf1DNB(stats::runif(1, -1, 1), stats::runif(1, 0, 5), stats::runif(1, 0, 5), stats::runif(1), stats::runif(1, 1, 2))))
})
