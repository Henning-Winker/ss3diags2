builds <- function(p) expect_silent(suppressWarnings(ggplot2::ggplot_build(p)))

test_that("observation plots build from FLR objects", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  ct <- ss3_catch(r); idx <- suppressMessages(ss3_indices(r)); fit <- suppressMessages(ss3_indices(r, what = "exp"))
  lfd <- ss3_lencomp(r); lfit <- ss3_lencomp(r, what = "exp")
  builds(ss3plot_data(catch = ct, indices = idx, lencomp = lfd))
  builds(ss3plot_catch(ct)); builds(ss3plot_catch(ct, type = "share")); builds(ss3plot_catch(ct, type = "line"))
  builds(ss3plot_indices(idx, fit = fit)); builds(ss3plot_indices(idx, combine = TRUE, log = TRUE))
  for(ty in c("pooled", "annual", "bubble", "mean")) builds(ss3plot_comps(lfd["Fishery_2"], fit = lfit["Fishery_2"],
                                                                         by = "unit", type = ty))
  builds(ss3plot_nsamp(lfd, by = "unit"))
  expect_error(ss3plot_comps(lfd, by = "sex"), "unit")
})

test_that("plots split by season and area", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  r$nseasons <- 2; r$nareas <- 2
  r$catch$Seas <- ifelse(r$catch$Yr %% 2 == 0, 1, 2); r$catch$Area <- ifelse(r$catch$Fleet %% 2 == 0, 2, 1)
  r$cpue$Seas <- ifelse(r$cpue$Fleet %% 2 == 0, 1, 2)
  r$lendbase$Seas <- ifelse(r$lendbase$Yr %% 2 == 0, 1, 2); r$lendbase$Area <- ifelse(r$lendbase$Fleet == 2, 2, 1)
  ct <- ss3_catch(r); lfd <- ss3_lencomp(r)
  expect_equal(dimnames(ct[[1]])$season, c("1", "2"))
  p <- ss3plot_catch(ct, by = c("season", "area")); builds(p)
  expect_equal(nlevels(p$data$group), 4)
  builds(ss3plot_indices(suppressMessages(ss3_indices(r)), by = "season"))
  builds(ss3plot_comps(lfd, by = c("unit", "area")))
  builds(ss3plot_data(catch = ct, lencomp = lfd, by = "season"))
})

test_that("stock plots build from FLStockR", {
  skip_no_data(); skip_if_not_installed("mvtnorm")
  e <- load_ss3("pac.hke"); r <- e$ss3phk
  stk <- ss3stk(r); stki <- ss3stk(r, output = "iters", iters = 200)
  builds(ss3plot_trj(stk)); builds(ss3plot_trj(stki)); builds(ss3plot_trj(list(MLE = stk, MVLN = stki)))
  builds(ss3plot_kobe(stki)); builds(ss3plot_kobe(stk))
  builds(ss3plot_production(stk)); builds(ss3plot_production(stk, type = "dyn", rel = TRUE))
})
