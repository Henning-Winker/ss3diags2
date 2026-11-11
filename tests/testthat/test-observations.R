test_that("ss3_fleets classifies fleets", {
  skip_no_data(); e <- load_ss3("natl.sma")
  fl <- ss3_fleets(e$ss3sma)
  expect_equal(nrow(fl), 18)
  expect_equal(sum(fl$type == "catch"), 12)
  expect_equal(fl[name == "CPUE_4"]$survey_units, "biomass")
})

test_that("ss3_catch matches ss3rep$catch", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  ct <- ss3_catch(r)
  expect_s4_class(ct, "FLQuants")
  expect_equal(length(ct), 12)
  raw <- r$catch[r$catch$Fleet == 1 & r$catch$Yr %in% r$startyr:r$endyr, ]
  expect_equal(as.numeric(ct[["Fishery_1"]]), raw$Obs[order(raw$Yr)])
  expect_equal(dimnames(ct[[1]])$year, as.character(r$startyr:r$endyr))
  kb <- ss3_catch(r, what = "kill_bio", long = TRUE)
  expect_true(all(kb$units == "t"))
})

test_that("ss3_indices builds FLIndexBiomass with SEs and timing", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  idx <- suppressMessages(ss3_indices(r))
  expect_s4_class(idx, "FLIndices")
  expect_equal(names(idx), paste0("CPUE_", 1:6))
  cp <- r$cpue[r$cpue$Fleet == 13, ]
  i1 <- idx[["CPUE_1"]]
  expect_equal(as.numeric(FLCore::index(i1))[!is.na(FLCore::index(i1))], cp$Obs[order(cp$Yr)])
  expect_equal(as.numeric(FLCore::index.var(i1))[!is.na(FLCore::index(i1))], cp$SE[order(cp$Yr)])
  expect_equal(unname(i1@range["startf"]), 0.5)
  expect_equal(FLCore::units(FLCore::index(idx[["CPUE_4"]])), "biomass")
  fit <- suppressMessages(ss3_indices(r, what = "exp"))
  expect_equal(as.numeric(FLCore::index(fit[[1]]))[!is.na(FLCore::index(fit[[1]]))], cp$Exp[order(cp$Yr)])
})

test_that("ss3_lencomp: sexes as unit, numbers = prop x nsamp", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  l <- ss3_lencomp(r)
  expect_equal(names(l), paste0("Fishery_", 1:5))
  expect_equal(dimnames(l[["Fishery_1"]])$unit, "unique")       # combined-sex samples
  expect_equal(dimnames(l[["Fishery_2"]])$unit, c("F", "M"))
  expect_equal(FLCore::quant(l[[1]]), "len")
  lb <- r$lendbase
  y <- lb$Yr[lb$Fleet == 2][1]
  s <- lb[lb$Fleet == 2 & lb$Yr == y, ]
  expect_equal(sum(l[["Fishery_2"]][, as.character(y)], na.rm = TRUE), sum(s$Obs*s$Nsamp_adj), tolerance = 1e-8)
  p <- ss3_lencomp(r, scale = "prop", sex = "combine")
  cs <- apply(p[["Fishery_2"]], 2, sum)
  expect_equal(unname(cs[!is.na(cs)]), rep(1, sum(!is.na(cs))), tolerance = 1e-6)
  lg <- ss3_lencomp(r, long = TRUE)
  expect_true(all(c("len", "year", "unit", "season", "area", "iter", "data", "obs", "exp", "nsamp") %in% names(lg)))
})

test_that("ss3_agecomp returns marginal ages by fleet", {
  skip_no_data(); e <- load_ss3("pac.hke"); r <- e$ss3phk
  a <- ss3_agecomp(r, scale = "prop")
  expect_equal(names(a), c("Fishery", "Acoustic_Survey"))
  expect_equal(FLCore::quant(a[[1]]), "age")
  cs <- apply(a[["Fishery"]], 2, sum)
  expect_equal(unname(cs[!is.na(cs)]), rep(1, sum(!is.na(cs))), tolerance = 1e-4)
  e2 <- load_ss3("gob.her")
  expect_equal(names(ss3_agecomp(e2$ss3her)), c("Fleet", "Acoustics", "Trapnet"))
})

test_that("ss2jabba3 returns annual catch FLQuants and FLIndices", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  j <- suppressMessages(ss2jabba3(r))
  expect_s4_class(j$catch, "FLQuants")
  expect_s4_class(j$indices, "FLIndices")
  tot <- tapply(r$catch$Obs, r$catch$Yr, sum)
  expect_equal(as.numeric(j$catch.total[, "2015"]), unname(tot["2015"]))
})

test_that("ss2flicc returns lfd on common bins and years plus catch by gear", {
  skip_no_data(); e <- load_ss3("natl.sma"); r <- e$ss3sma
  f <- ss2flicc(r)
  expect_equal(names(f$lfd), paste0("Fishery_", 1:5))
  expect_true(all(sapply(f$lfd, function(q) identical(dimnames(q)$len, dimnames(f$lfd[[1]])$len))))
  expect_true(all(sapply(f$lfd, function(q) identical(dimnames(q)$year, as.character(f$years)))))
  expect_false(anyNA(unlist(lapply(f$lfd, as.numeric))))
  expect_equal(names(f$catch_by_gear), names(f$lfd))
  fa <- ss2flicc(r, years = "all")
  expect_true(anyNA(as.numeric(fa$lfd[["Fishery_5"]])))
})
