test_that("ss3mvln matches ss3diags::SSdeltaMVLN at the MLE and in distribution", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  for(nm in c("natl.sma", "pac.hke", "gob.her")){
    e <- load_ss3(nm); r <- get(ls(e, pattern = "^ss3")[1], e)
    a <- ss3mvln(r, mc = 2000, method = "annual")
    o <- ss3diags::SSdeltaMVLN(r, mc = 2000, plot = FALSE, verbose = FALSE)
    m <- merge(a$mle, o$mle, by = "year")
    expect_equal(a$method, "annual")
    for(v in c("stock", "harvest", "F", "SSB", "Recr", "Catch"))
      expect_equal(m[[paste0(v, ".x")]], m[[paste0(v, ".y")]], tolerance = 1e-10, label = paste(nm, v))
    expect_equal(as.numeric(a$refpts), o$refpts$value, tolerance = 1e-10)
    # draws: same marginals; correlations averaged over all neighbouring years (SSdeltaMVLN recycles labels)
    ly <- max(a$kb$year)
    expect_equal(stats::median(a$kb[year == ly]$stock), stats::median(o$kb$stock[o$kb$year == ly]), tolerance = 0.01)
    expect_equal(stats::median(a$kb[year == ly]$harvest), stats::median(o$kb$harvest[o$kb$year == ly]), tolerance = 0.01)
  }
})

test_that("ss3stk builds MLE and MVLN FLStockR with refpts", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  e <- load_ss3("natl.sma"); r <- e$ss3sma
  stk <- ss3stk(r, blim = 0.4)
  expect_s4_class(stk, "FLStockR")
  hat <- r$derived_quants
  expect_equal(as.numeric(FLCore::ssb(stk)[, "2015"]), hat$Value[hat$Label == "SSB_2015"])
  expect_equal(as.numeric(FLCore::fbar(stk)[, "2015"]), ss3mvln(r, mc = 2)$mle[year == 2015]$F)
  expect_true(all(c("Ftgt", "Btgt", "MSY", "B0", "R0", "Blim") %in% dimnames(FLCore::refpts(stk))$params))
  expect_s4_class(attr(stk, "se.log"), "FLQuants")
  stki <- ss3stk(r, output = "iters", iters = 50)
  expect_equal(dim(FLCore::ssb(stki))[6], 50)
  expect_equal(dim(FLCore::refpts(stki))[2], 50)   # joint MVLN: refpts with iterations
  expect_true(nrow(ss3_eqcurve(r)) > 10)
})

test_that("ss3mvln joint: MLE as annual, coherent draws, CoVar correlations, refpts with iterations", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  e <- load_ss3("natl.sma"); r <- e$ss3sma
  j <- ss3mvln(r, mc = 2000); a <- ss3mvln(r, mc = 2, method = "annual")
  expect_s3_class(j, "ss3mvln"); expect_equal(j$method, "joint")
  m <- merge(j$mle, a$mle, by = "year")
  for(v in c("stock", "harvest", "F", "SSB", "Recr", "Catch"))
    expect_equal(m[[paste0(v, ".x")]], m[[paste0(v, ".y")]], tolerance = 1e-5, label = v)
  expect_equal(as.numeric(j$refpts), as.numeric(a$refpts), tolerance = 1e-8)
  expect_equal(dim(j$refpts.iter)[2], 2000)
  ly <- r$endyr; k <- j$kb[year == ly]
  cv <- r$CoVar; rho <- cv$corr[cv$label.i == paste0("Bratio_", ly) & cv$label.j == paste0("F_", ly)]
  expect_equal(stats::cor(log(k$stock), log(k$harvest)), rho, tolerance = 0.05)
  expect_equal(stats::median(k$stock), j$mle[year == ly]$stock, tolerance = 0.01)
  w <- data.table::dcast(j$kb[year %in% c(ly - 1, ly)], iter ~ year, value.var = "SSB")
  expect_gt(stats::cor(w[[2]], w[[3]]), 0.9)                          # iterations are coherent trajectories
  expect_error(ss3mvln(load_ss3("natl.sma")$ss3sma[setdiff(names(r), "CoVar")]), "CoVar")
})

test_that("ss3stk: joint iterations with iteration refpts, rel, thin; ss3_flqs", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  e <- load_ss3("pac.hke"); r <- e$ss3phk
  s <- ss3stk(r, output = "iters", iters = 200, thin = 2)
  expect_equal(dim(FLCore::ssb(s))[6], 100)
  expect_equal(dim(FLCore::refpts(s))[2], 100)
  q <- ss3_flqs(s)
  expect_s4_class(q, "FLQuants"); expect_equal(names(q), c("SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt"))
  mv <- ss3mvln(r, mc = 200)
  expect_equal(as.numeric(q[["B/Btgt"]][, "2019", , , , 1]), mv$kb[year == 2019 & iter == 1]$stock, tolerance = 1e-8)
  sr <- ss3stk(mv, rel = TRUE)
  expect_equal(as.numeric(FLCore::refpts(sr)["Btgt", 1]), 1)
  expect_equal(as.numeric(FLCore::ssb(sr)[, "2019", , , , 1]), mv$kb[year == 2019 & iter == 1]$stock, tolerance = 1e-8)
  sa <- ss3stk(ss3mvln(r, mc = 50, method = "annual"))
  expect_s4_class(attr(sa, "status"), "FLQuants")
  expect_silent(ggplot2::ggplot_build(ss3plot_kobe(s)))
})
