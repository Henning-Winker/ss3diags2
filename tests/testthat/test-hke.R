# European hake GSAs 17-18: two sexes, conditional age-at-length
hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))

test_that("conditional age-at-length residuals, Francis weights and runs on real data", {
  a <- ss3_resid(hke, "con", long = TRUE)
  i <- ss3_resid(hke, "con", long = TRUE, nsamp = "input")
  expect_gt(sum(is.finite(i$data)), sum(is.finite(a$data)))
  expect_setequal(unique(i$unit), c("F", "M"))
  lh <- suppressMessages(ss3_lhpar(hke))
  expect_silent(ggplot2::ggplot_build(ss3plot_caal(i, fleet = "ITA_LLS_18", years = 2020:2021, lhpar = lh)))
  f <- ss3_francis(hke, "con")
  expect_equal(f$fleet, c("ITA_OTB_18", "ITA_LLS_18"))
  expect_equal(f$var_adj, c(0.01, 0.11))
  skip_if_not_installed("r4ss")
  b <- vapply(5:6, function(k) unname(r4ss::SSMethod.Cond.TA1.8(hke, fleet = k, plotit = FALSE)[1]), 1)
  expect_equal(f$w, b, tolerance = 1e-8)
})

test_that("impact by sex reproduces SS3 SSB and dynamic B0 (sexes with different M)", {
  ck <- ss3_impact(hke)$check
  expect_equal(ck$projected, ck$ss3, tolerance = 1e-8)
  expect_equal(ck$dynB0, ck$ss3_dynB0, tolerance = 1e-4)
})
