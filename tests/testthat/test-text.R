test_that("text = FALSE drops title, subtitle and caption but keeps legends", {
  data(ss3simple, package = "ss3diags2", envir = environment())
  r <- ss3_resid(ss3simple)
  p <- ss3plot_resid(r, lhpar = ss3_lhpar(ss3simple), text = FALSE)
  expect_null(p$labels$title); expect_null(p$labels$subtitle); expect_null(p$labels$caption)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  g <- ss3plot_growth(ss3simple, lfd = ss3_lencomp(ss3simple), text = FALSE)
  expect_true(any(vapply(g$layers, function(l) inherits(l$geom, "GeomText"), TRUE)))   # in-plot labels stay
  expect_null(ss3plot_kobe(ss3stk(ss3simple, output = "iters", iters = 50), text = FALSE)$labels$caption)
  expect_identical(ss3plot_sr(ss3simple, text = TRUE)$labels$title, ss3plot_sr(ss3simple)$labels$title)
})
