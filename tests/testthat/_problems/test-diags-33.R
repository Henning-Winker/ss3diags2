# Extracted from test-diags.R:33

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ss3diags2", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(retro.sma, package = "ss3diags2", envir = environment())
e <- ss3_ensemble(list(a = retro.sma$replist1, b = retro.sma$replist2), weights = c(3, 1), iters = 40)
