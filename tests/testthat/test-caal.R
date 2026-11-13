test_that("OSA residuals of conditional age-at-length data are N(0, 1) under the model", {
  data(ss3simple, package = "ss3diags2", envir = environment())
  r <- fake_caal(ss3simple)
  d <- ss3_resid(r, quant = "con", long = TRUE)
  expect_true(all(c("lbin", "age") %in% names(d)))
  x <- d[is.finite(data)]$data
  expect_lt(abs(mean(x)), 0.1); expect_lt(abs(sd(x) - 1), 0.1)
  q <- ss3_resid(r, quant = "con")
  expect_s4_class(q, "FLQuants"); expect_equal(names(q), "SURVEY1")
  expect_silent(ggplot2::ggplot_build(ss3plot_caal(d, lhpar = ss3_lhpar(r))))
  expect_silent(ggplot2::ggplot_build(ss3plot_resid(q)))
})

test_that("runs tests, joint residuals and hindcast work on conditional age-at-length", {
  data(ss3simple, package = "ss3diags2", envir = environment())
  r <- fake_caal(ss3simple)
  ru <- ss3_runs(r, quant = "con")
  expect_equal(ru$index, "SURVEY1"); expect_true(ru$test %in% c("Passed", "Failed"))
  expect_s3_class(ss3plot_runs(r, quant = "con"), "ggplot")
  expect_s3_class(ss3plot_jres(r, quant = "con"), "ggplot")
  # hindcast: peels drop the last years into the ghost table
  peel <- function(k){ x <- r; cb <- x$condbase; ey <- r$endyr - k
    x$ghostcondbase <- cb[cb$Yr > ey, ]; x$ghostcondbase$Fleet <- -x$ghostcondbase$Fleet
    x$condbase <- cb[cb$Yr <= ey, ]; x$Retro_year <- ey; x }
  hc <- ss3_hindcast(c(list(r), lapply(1:4, peel)), quant = "con")
  expect_equal(hc$mase$index[1], "SURVEY1")
  expect_true(is.finite(hc$mase$MASE[1]))
})
