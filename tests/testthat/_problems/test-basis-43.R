# Extracted from test-basis.R:43

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ss3diags2", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
data(ss3simple, package = "ss3diags2")
b <- function(r, Fref = NULL) ss3diags2:::.ss3_basis(r, Fref)
setq <- function(r, l, v){ r$derived_quants$Value[r$derived_quants$Label == l] <- v; r }
gq <- function(r, l) r$derived_quants$Value[r$derived_quants$Label == l]
ey <- ss3simple$endyr
rebratio <- function(r, den){
  for(l in grep("^Bratio_", r$derived_quants$Label, value = TRUE))
    r <- setq(r, l, gq(r, sub("Bratio", "SSB", l))/den)
  r
}
msy_as <- function(r, what){                     # MSY option set to F(SPR) or F(Btgt), Bratio = B/Bmsy
  r <- setq(r, "SSB_MSY", gq(r, paste0("SSB_", what))); r <- setq(r, "annF_MSY", gq(r, paste0("annF_", what)))
  r <- setq(r, "Dead_Catch_MSY", gq(r, paste0("Dead_Catch_", what)))
  rebratio(r, gq(r, "SSB_MSY"))
}
stat <- function(r, Fref = NULL) unlist(ss3mvln(r, mc = 0, Fref = Fref)$mle[year == ey, .(stock, harvest)])

# test -------------------------------------------------------------------------
skip_if_not_installed("ss3diags")
for(r in list(ss3simple, msy_as(ss3simple, "SPR"), msy_as(ss3simple, "Btgt"),
                rebratio(ss3simple, gq(ss3simple, "SSB_MSY")))){
    for(Fref in list(NULL, "SPR", "Btgt")){
      d <- ss3diags::SSdeltaMVLN(r, mc = 10, plot = FALSE, verbose = FALSE, Fref = Fref)$mle
      expect_equal(stat(r, Fref), unlist(d[d$year == ey, c("stock", "harvest")]), tolerance = 1e-6, ignore_attr = TRUE)
    }
  }
