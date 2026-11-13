test_that("Francis TA1.8 weights match r4ss::SSMethod.TA1.8", {
  data(ss3simple, package = "ss3diags2", envir = environment())
  w <- ss3_francis(ss3simple)
  expect_equal(w$fleet, c("FISHERY", "SURVEY1"))
  expect_equal(w$var_adj, c(1, 1))
  expect_true(all(w$lo < w$w & w$w < w$hi))
  skip_if_not_installed("r4ss")
  for(q in c("len", "age")){
    a <- ss3_francis(ss3simple, q)
    b <- vapply(1:2, function(f) unname(r4ss::SSMethod.TA1.8(ss3simple, type = q, fleet = f, plotit = FALSE)[1]), 1)
    expect_equal(a$w, b, tolerance = 1e-8)
  }
  expect_s3_class(ss3plot_francis(w), "ggplot")
  expect_silent(ggplot2::ggplot_build(ss3plot_francis(ss3_francis(fake_caal(ss3simple), "con"))))
})
