test_that("catch options, probabilities and plot from forecast scenarios", {
  fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
  fs <- ss3_fwd(fw, output = "iters", iters = 300, Fref = "Btgt", blim = 0.25, bpa = 0.5)
  co <- ss3_catchopts(fs, bases = c("F0", "Fmsy", "Fsq", "P(Blim)"), btrigger = 0.5)
  expect_equal(attr(co, "year"), 2027)
  expect_equal(co$basis, c("F0", "Fmsy", "Fsq", "P(Blim)"))
  cc <- co[["catch 2027"]]
  expect_true(cc[1] < 100)                                   # F0 scenario
  expect_true(cc[2] > cc[1])
  sm <- attr(co, "scenarios")
  # Fmsy option: F/Fmsy = 1 (SSB above Btrigger, no reduction) unless outside the scenarios
  if(is.finite(cc[2])) expect_equal(co[["F2027/Fmsy"]][2], 1, tolerance = 0.01)
  expect_s3_class(ss3plot_catchopts(co), "ggplot")
  pr <- ss3_fwdprobs(fs, btrigger = 0.5)
  expect_named(pr, c("F<Ftgt", "SSB>Btgt", "green", "SSB>Btrigger", "SSB>Blim"))
  expect_equal(pr$`SSB>Blim`$scenario, names(fs))
  expect_true(all(unlist(pr$green[, -1]) >= 0 & unlist(pr$green[, -1]) <= 100))
  # risk from the ratio draws equals the share of B/Btgt draws below Blim/Btgt
  s <- fs[[1]]; r <- as.numeric(attr(s, "status")$stock[, "2028"])
  expect_equal(ss3diags2:::.p_below(s, 2028, "Blim"), mean(r < 0.25))
})
