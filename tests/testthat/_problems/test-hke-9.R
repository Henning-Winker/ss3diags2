# Extracted from test-hke.R:9

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ss3diags2", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))

# test -------------------------------------------------------------------------
a <- ss3_resid(hke, "con", long = TRUE)
i <- ss3_resid(hke, "con", long = TRUE, nsamp = "input")
expect_gt(sum(is.finite(i$data)), sum(is.finite(a$data)))
expect_setequal(unique(i$unit), c("F", "M"))
expect_silent(ggplot2::ggplot_build(ss3plot_caal(i, fleet = "ITA_LLS_18", years = 2020:2021, lhpar = ss3_lhpar(hke))))
