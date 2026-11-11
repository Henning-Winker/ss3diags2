test_that("ss3_mvln matches ss3diags::SSdeltaMVLN at the MLE and in distribution", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  for(nm in c("natl.sma", "pac.hke", "gob.her")){
    e <- load_ss3(nm); r <- get(ls(e, pattern = "^ss3")[1], e)
    a <- ss3_mvln(r, mc = 2000)
    o <- ss3diags::SSdeltaMVLN(r, mc = 2000, plot = FALSE, verbose = FALSE)
    m <- merge(a$mle, o$mle, by = "year")
    for(v in c("stock", "harvest", "F", "SSB", "Recr", "Catch"))
      expect_equal(m[[paste0(v, ".x")]], m[[paste0(v, ".y")]], tolerance = 1e-10, label = paste(nm, v))
    expect_equal(as.numeric(a$refpts), o$refpts$value, tolerance = 1e-10)
    # draws: same marginals; correlations averaged over all neighbouring years (SSdeltaMVLN recycles labels)
    ly <- max(a$kb$year)
    expect_equal(stats::median(a$kb[year == ly]$stock), stats::median(o$kb$stock[o$kb$year == ly]), tolerance = 0.01)
    expect_equal(stats::median(a$kb[year == ly]$harvest), stats::median(o$kb$harvest[o$kb$year == ly]), tolerance = 0.01)
  }
})

test_that("ss3_FLStockR builds MLE and MVLN FLStockR with refpts", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  e <- load_ss3("natl.sma"); r <- e$ss3sma
  stk <- ss3_FLStockR(r, blim = 0.4)
  expect_s4_class(stk, "FLStockR")
  hat <- r$derived_quants
  expect_equal(as.numeric(FLCore::ssb(stk)[, "2015"]), hat$Value[hat$Label == "SSB_2015"])
  expect_equal(as.numeric(FLCore::fbar(stk)[, "2015"]), ss3_mvln(r, mc = 2)$mle[year == 2015]$F)
  expect_true(all(c("Ftgt", "Btgt", "MSY", "B0", "R0", "Blim") %in% dimnames(FLCore::refpts(stk))$params))
  expect_s4_class(attr(stk, "se.log"), "FLQuants")
  stki <- ss3_FLStockR(r, iters = 50)
  expect_equal(dim(FLCore::ssb(stki))[6], 50)
  expect_s4_class(attr(stki, "status"), "FLQuants")
  expect_true(nrow(ss3_eqcurve(r)) > 10)
})
