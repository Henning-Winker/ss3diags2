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

test_that("Kobe II strategy matrices (ICCAT, IOTC)", {
  fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
  fs <- ss3_fwd(fw, output = "iters", iters = 200, Fref = "Btgt", blim = 0.25, bpa = 0.5)
  k1 <- ss3_k2sm(fs, "iccat")
  expect_named(k1, c("green", "SSB>=SSBtgt", "F<=Ftgt"))
  g <- as.matrix(k1$green[, -(1:2)]); b <- as.matrix(k1$`SSB>=SSBtgt`[, -(1:2)])
  expect_true(all(g <= b + 1))                               # green is a subset of SSB >= SSBtgt
  k2 <- ss3_k2sm(fs, "iotc", flim = 1.4)
  expect_true(all(c("SSB<SSBtgt", "F>Ftgt", "SSB<SSBlim", "F>Flim") %in% names(k2)))
  expect_equal(as.matrix(k2$`SSB<SSBtgt`[, -(1:2)]), 100 - b, tolerance = 1)
  expect_s3_class(ss3plot_k2sm(k1), "ggplot")
})

test_that("ICES and GFCM risk matrices", {
  fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
  fs <- ss3_fwd(fw, output = "iters", iters = 200, Fref = "Btgt", blim = 0.25, bpa = 0.5)
  k3 <- ss3_k2sm(fs, "ices", btrigger = 0.5)
  expect_named(k3, c("SSB<Blim", "SSB<Btrigger", "F>Fmsy", "green"))
  k4 <- ss3_k2sm(fs, "gfcm")
  expect_named(k4, c("SSB<Blim", "SSB<Bpa", "F>Ftgt", "green"))
  st <- attr(k4, "stars")
  expect_true(all(c("C", "F/Ftgt", "B/Btgt", "B/Bpa", "B/Blim") %in% names(st)))
  expect_equal(st$`B/Blim`, round(st$`B/Btgt`/0.25, 2), tolerance = 0.02)
})
