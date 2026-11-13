# Extracted from test-caal.R:10

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ss3diags2", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(ss3simple, package = "ss3diags2", envir = environment())
r <- fake_caal(ss3simple)
d <- ss3_resid(r, quant = "con", long = TRUE)
expect_true(all(c("lbin", "age") %in% names(d)))
x <- d[is.finite(data)]$data
expect_lt(abs(mean(x)), 0.1)
expect_lt(abs(sd(x) - 1), 0.1)
q <- ss3_resid(r, quant = "con")
expect_s4_class(q, "FLQuants")
expect_equal(names(q), "SURVEY1")
expect_silent(ggplot2::ggplot_build(ss3plot_caal(d, lhpar = ss3_lhpar(r))))
