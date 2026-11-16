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

test_that("sex ratio plot and CAAL means", {
  lf <- ss3_lencomp(hke); lx <- ss3_lencomp(hke, what = "exp")
  p <- ss3plot_sexratio(lf, fit = lx)
  expect_true(all(p$data$pF >= 0 & p$data$pF <= 1))
  expect_false("MNE_OTB_18" %in% p$data$name)
  expect_silent(ggplot2::ggplot_build(ss3plot_sexratio(lf, fit = lx, type = "year")))
  lh <- suppressMessages(ss3_lhpar(hke))
  expect_equal(c(lh["richards", "F"]), hke$parameters["Richards_Fem_GP_1", "Value"])
  g <- ss3diags2:::.growth_curves(lh, maxage = 12)
  eg <- hke$endgrowth[hke$endgrowth$Seas == 1 & hke$endgrowth$Sex == 1 & hke$endgrowth$int_Age == 5, "Len_Mid"][1]
  expect_equal(stats::approx(g[unit == "F"]$age, g[unit == "F"]$L, 5.5)$y, eg, tolerance = 1e-3)
})
