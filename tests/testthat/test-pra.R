# season-, sex- and area-structured model: northern shrimp pra.27.3a4a
data(ss3pra, package = "ss3diags2")
g3 <- list(Denmark = 1:2, Norway = 3:4, Sweden = 5:6)

test_that("observations keep season, unit and area structure", {
  ct <- ss3_catch(ss3pra)
  expect_equal(dim(ct[[1]])[4], 4)                                  # seasons
  expect_equal(dimnames(ct[["DK4a"]])$area, "2")                   # fleet area
  lf <- ss3_lencomp(ss3pra)
  expect_true(all(c("F", "M") %in% dimnames(lf[["NORSURVEY3a"]])$unit))
  expect_true(any(grepl("discard", names(lf))))
})

test_that("ss3_fa returns F-at-age by fleet with age, unit, season and area", {
  fa <- ss3_fa(ss3pra)
  expect_equal(names(fa), paste0(rep(c("DK", "NOR", "SWE"), each = 2), c("3a", "4a")))
  expect_equal(dim(fa[[1]])[c(1, 3, 4)], c(9, 2, 4))
  fan <- ss3_fa(ss3pra, annual = TRUE)
  expect_equal(c(fan[[1]][, "2020", 1]), unname(apply(fa[[1]][, "2020", 1], 1, sum)*0.25), tolerance = 1e-10)
})

test_that("partial F by catch adds up to the reported F", {
  pf <- ss3_partialF(ss3pra)
  expect_equal(c(Reduce(`+`, lapply(pf, c))), c(attr(pf, "total")), tolerance = 1e-8)
  pg <- ss3_partialF(ss3pra, groups = g3)
  expect_equal(names(pg), names(g3))
  fa <- ss3_partialF(ss3pra, groups = "area", method = "fbar")
  expect_equal(names(fa), c("area 1", "area 2"))
  expect_s3_class(ss3plot_partialF(pf), "ggplot")
  expect_s3_class(ss3plot_partialF(ss3_partialF(ss3pra, method = "apical", season = TRUE)), "ggplot")
})

test_that("partial impact reproduces SS3 SSB and dynamic B0", {
  imp <- ss3_impact(ss3pra, groups = g3)
  ck <- imp$check
  expect_equal(ck$dynB0, ck$ss3_dynB0, tolerance = 1e-4)          # no fishing = SS3 dynamic B0
  late <- ck$year >= 2000
  expect_equal(ck$projected[late], ck$ss3[late], tolerance = 1e-4)  # all fleets = SS3 SSB
  im <- imp$impact[year == 2024]
  b0 <- imp$ssb[scenario == "no fishing" & year == 2024]$ssb
  expect_equal(sum(im$impact), b0 - ck[year == 2024]$projected, tolerance = 1e-8)
  m <- ss3_impact(ss3pra, groups = g3, init = "model")
  expect_equal(m$check$projected, m$check$ss3, tolerance = 1e-8)
  expect_s3_class(ss3plot_impact(imp), "ggplot")
  expect_s3_class(ss3plot_impact(ss3_impact(ss3pra, type = "remove"), type = "line"), "ggplot")
})

test_that("selectivity, stock by area and stock-recruitment", {
  s <- ss3_sel(ss3pra)
  expect_equal(length(s), 8)
  expect_true(all(c("F", "M") %in% dimnames(s[[1]])$unit))
  ba <- ss3_byarea(ss3pra)
  expect_equal(dimnames(ba$SSB)$area, c("1", "2"))
  tot <- data.table::as.data.table(ss3pra$timeseries)[Yr == 2020 & Seas == 1, sum(SpawnBio)]
  expect_equal(sum(ba$SSB[, "2020"]), tot)
  expect_s3_class(ss3plot_sel(s, ret = ss3_sel(ss3pra, "ret")), "ggplot")
  expect_s3_class(ss3plot_byarea(ba), "ggplot")
  expect_s3_class(ss3plot_sr(ss3pra), "ggplot")
  expect_s3_class(ss3plot_fa(ss3_fa(ss3pra)), "ggplot")
  expect_s3_class(ss3plot_catch(ss3_catch(ss3pra), by = c("season", "area")), "ggplot")
})

test_that("plots select fleets and units", {
  rs <- ss3_resid(ss3pra)
  lh <- ss3_lhpar(ss3pra)
  p <- ss3plot_resid(rs, lhpar = lh, fleets = "DK3a")
  expect_setequal(unique(as.character(p$data$name)), c("DK3a.discard", "DK3a.retained"))
  p <- ss3plot_resid(rs, units = "F")
  expect_true(all(p$data$unit == "F"))
  expect_true(all(grepl("SURVEY|SWE3a", unique(as.character(p$data$name)))))
  lf <- ss3_lencomp(ss3pra); lx <- ss3_lencomp(ss3pra, what = "exp")
  expect_s3_class(ss3plot_comps(lf, fit = lx, fleets = c("NORSURVEY3a", "NORSURVEY4a"), units = "M"), "ggplot")
  expect_s3_class(ss3plot_ridges(lf, fit = lx, type = "both", fleets = "NORSURVEY3a", units = "F", lhpar = lh), "ggplot")
  expect_error(ss3plot_resid(rs, fleets = "nope"), "no data")
})
