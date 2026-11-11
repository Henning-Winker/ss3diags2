test_that("anglerfish: partitions, sex-specific surveys and common length bins", {
  data(ss3ank, package = "ss3diags2", envir = environment())
  l <- ss3_lencomp(ss3ank)
  expect_equal(names(l), c("FL1_Trawls.discard", "FL1_Trawls.retained", "FL2_Gillnets.discard",
                           "FL2_Gillnets.retained", "FR_IE_IBTS", "IE_IAMS"))
  expect_equal(dimnames(l[["FR_IE_IBTS"]])$unit, c("F", "M"))
  expect_equal(dimnames(l[["FL1_Trawls.retained"]])$unit, "unique")
  expect_true(all(sapply(l, function(q) identical(dimnames(q)$len, as.character(ss3ank$lbins)))))
  # sampled years complete (0 for compressed tail bins), unsampled years NA
  q <- l[["FL1_Trawls.discard"]]; yy <- colSums(q[, , 1, 1, 1, 1, drop = TRUE], na.rm = FALSE)
  expect_false(anyNA(q[, names(yy)[is.finite(yy)]]))
  lb <- ss3ank$lendbase; s <- lb[lb$Fleet == 3 & lb$Yr == 2010, ]
  expect_equal(sum(l[["FR_IE_IBTS"]][, "2010"]), sum(s$Obs*s$Nsamp_adj), tolerance = 1e-8)
})

test_that("anglerfish: input SE, discards, newer r4ss catch columns", {
  data(ss3ank, package = "ss3diags2", envir = environment())
  idx <- ss3_indices(ss3ank)
  cp <- ss3ank$cpue[ss3ank$cpue$Fleet == 3, ]
  v <- as.numeric(FLCore::index.var(idx[["FR_IE_IBTS"]])); v <- v[!is.na(v)]
  expect_equal(v, cp$SE_input[order(cp$Yr)])
  d <- ss3_discards(ss3ank)
  expect_equal(names(d), c("FL1_Trawls", "FL2_Gillnets"))
  expect_equal(as.numeric(d[[1]][, "2004"]), 676)
  kb <- ss3_catch(ss3ank, what = "kill_bio", long = TRUE)          # dead_bio in this r4ss version
  expect_true(all(kb$data >= ss3_catch(ss3ank, what = "ret_bio", long = TRUE)$data - 1e-6))
})

test_that("anglerfish: FLStockR from the time series when derived quantities are 0", {
  data(ss3ank, package = "ss3diags2", envir = environment())
  stk <- suppressMessages(ss3stk(ss3ank, btgt = 0.4))
  ts <- ss3ank$timeseries
  expect_equal(as.numeric(FLCore::ssb(stk)[, "2025"]), ts$SpawnBio[ts$Era == "TIME" & ts$Yr == 2025])
  expect_equal(as.numeric(FLCore::refpts(stk)["B0"]), ts$SpawnBio[ts$Era == "VIRG"])
  expect_equal(as.numeric(FLCore::refpts(stk)["Btgt"]), 0.4*ts$SpawnBio[ts$Era == "VIRG"])
  expect_error(ss3stk(ss3ank, output = "iters", iters = 10), "CoVar")
  expect_silent(suppressMessages(ggplot2::ggplot_build(ss3plot_trj(stk))))
})
