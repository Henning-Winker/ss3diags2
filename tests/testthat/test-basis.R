# reference point basis and labels: logic of ss3diags::SSdeltaMVLN, labels named by benchmark
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

test_that("labels name the benchmark behind the targets", {
  x <- b(ss3simple); expect_equal(c(x$Bref, x$Fref, x$Cref), c("B40", "Fmsy", "MSY"))
  x <- b(rebratio(ss3simple, gq(ss3simple, "SSB_MSY"))); expect_equal(c(x$Bref, x$Fref, x$Cref), c("Bmsy", "Fmsy", "MSY"))
  x <- b(msy_as(ss3simple, "SPR")); expect_equal(c(x$Bref, x$Fref, x$Cref), c("Bspr50", "Fspr50", "Yspr50"))
  x <- b(msy_as(ss3simple, "Btgt")); expect_equal(c(x$Bref, x$Fref, x$Cref), c("B40", "F40", "Y40"))
  pra <- readRDS(system.file("extdata", "pra.27.3a4a", "pra.27.3a4a.rds", package = "ss3diags2"))
  x <- b(pra); expect_equal(c(x$Bref, x$Fref, x$Cref), c("B30", "F30", "Y30"))
})

test_that("F report bases as ratios are recognised", {
  fy <- grep("^F_[0-9]", ss3simple$derived_quants$Label, value = TRUE)
  fr <- function(den, basis){ r <- ss3simple; for(l in fy) r <- setq(r, l, gq(r, l)/den); r$F_std_basis <- basis; r }
  d <- fr(gq(ss3simple, "annF_Btgt"), "(F)/(F_at_B40%);_with_F=sum(full_Fs)")
  expect_equal(b(d)$Fquant, "Btgt"); expect_equal(stat(d), stat(ss3simple, "Btgt"), tolerance = 1e-8)
  s <- fr(gq(ss3simple, "annF_SPR"), "(F)/(F50%SPR);_with_F=sum(full_Fs)")
  expect_equal(b(s)$Fref, "Fspr50"); expect_equal(stat(s), stat(ss3simple, "SPR"), tolerance = 1e-8)
})

test_that("stock status matches ss3diags::SSdeltaMVLN", {
  skip_if_not_installed("ss3diags")
  for(r in list(ss3simple, msy_as(ss3simple, "SPR"), msy_as(ss3simple, "Btgt"),
                rebratio(ss3simple, gq(ss3simple, "SSB_MSY")))){
    for(Fref in list(NULL, "SPR", "Btgt")){
      d <- ss3diags::SSdeltaMVLN(r, mc = 10, plot = FALSE, verbose = FALSE, Fref = Fref)$mle
      expect_equal(stat(r, Fref), unlist(d[d$year == ey, c("stock", "harvest")]), tolerance = 1e-6, ignore_attr = TRUE)
    }
  }
})
