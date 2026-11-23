test_that("forecast scenarios from SSsummarize output: FLStocks, intermediate years, plot and table", {
  fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
  hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
  fs <- ss3_fwd(fw, Fref = "Btgt", bpa = 0.5, blim = 0.25)
  expect_equal(names(fs), c("Fadv", "Fsq", "F0"))
  expect_equal(ss3diags2:::.fwd_nint(fs, 2024), 2)
  expect_equal(unname(round(fs[[1]]@refpts["Blim"][[1]])), 23807)
  stk <- ss3stk(hke, Fref = "Btgt", bpa = 0.5, blim = 0.25)
  expect_equal(c(FLCore::ssb(stk)[, "2024"]), c(FLCore::ssb(fs$Fadv)[, "2024"]))      # Fadv run = reference model
  expect_s3_class(ss3plot_fwd(fs, stk = stk, ci = FALSE), "ggplot")
  tb <- ss3_fwdtable(fs)
  expect_equal(names(tb)[2:4], c("C_2027", "F_2027", "SSB_2028"))
  expect_true(tb[Basis == "F0"]$SSB_2028 > tb[Basis == "Fadv"]$SSB_2028)
  expect_true(all(tb[[6]] >= 0 & tb[[6]] <= 1))
})

test_that("ensemble forecast: weights, resampling beyond the draws, per-draw Blim", {
  hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
  fa <- ss3_fwd(list(Fadv = hke), output = "iters", iters = 60, Fref = "Btgt", bpa = 0.5, blim = 0.25)
  fb <- ss3_fwd(list(Fadv = hke), output = "iters", iters = 60, Fref = "Btgt", bpa = 0.5, blim = 0.20)
  fe <- ss3_fwdensemble(list(A = fa, B = fb), weights = c(0.6, 0.4), iters = 100)   # 120 > 60 draws of A
  expect_equal(attr(fe, "weights")$iters, c(120, 80))
  expect_equal(FLCore::dims(fe$Fadv)$iter, 200)
  bl <- as.numeric(FLCore::refpts(fe$Fadv)["Blim", ])
  expect_length(bl, 200)
  expect_true(all(is.finite(bl)))
  tb <- ss3_fwdtable(fe)
  y1 <- grep("^SSB_", names(tb), value = TRUE)[1]
  p <- mean(as.numeric(FLCore::ssb(fe$Fadv)[, sub("SSB_", "", y1)]) < bl)
  expect_equal(tb[[6]], round(p, 3))
  expect_s3_class(ss3plot_fwd(fe), "ggplot")
})
