test_that("OSA residuals are standard normal for simulated multinomial and DM compositions", {
  set.seed(3); p <- stats::dnorm(1:30, 15, 4); p <- p/sum(p)
  rm <- unlist(replicate(200, ss3diags2:::.osa_comp(as.numeric(stats::rmultinom(1, 200, p)), p), simplify = FALSE))
  expect_lt(abs(mean(rm, na.rm = TRUE)), 0.05); expect_lt(abs(stats::sd(rm, na.rm = TRUE) - 1), 0.05)
  a0 <- 300
  rd <- unlist(replicate(200, { pp <- stats::rgamma(30, a0*p); pp <- pp/sum(pp)
    ss3diags2:::.osa_comp(as.numeric(stats::rmultinom(1, 200, pp)), p, a0) }, simplify = FALSE))
  expect_lt(abs(mean(rd, na.rm = TRUE)), 0.05); expect_lt(abs(stats::sd(rd, na.rm = TRUE) - 1), 0.05)
})

test_that("ss3_resid uses the DM precision implied by SS3 and returns FLQuants by fleet", {
  data(ss3ank, package = "ss3diags2", envir = environment())
  d <- ss3_resid(ss3ank, long = TRUE)
  lb <- ss3ank$lendbase; s <- lb[lb$Fleet == 1 & lb$Yr == 2003 & lb$Part == 1, ][1, ]
  th <- exp(ss3ank$parameters["ln(DM_theta)_Len_P1", "Value"])
  expect_equal(d[partname == "FL1_Trawls.discard" & year == 2003]$a0[1], th*s$Nsamp_adj, tolerance = 0.01)
  r <- ss3_resid(ss3ank)
  expect_s4_class(r, "FLQuants")
  expect_equal(attr(r, "type"), "osa")
  expect_equal(dimnames(r[["FR_IE_IBTS"]])$unit, c("F", "M"))
  pe <- ss3_resid(ss3ank, type = "pearson", long = TRUE)
  expect_equal(sum(is.finite(pe$data)), sum(is.finite(lb$Pearson)))
})

test_that("ss3_lhpar returns von Bertalanffy parameters by unit and residual plots build", {
  data(ss3ank, package = "ss3diags2", envir = environment())
  lh <- ss3_lhpar(ss3ank)
  expect_equal(dimnames(lh)$unit, c("F", "M"))
  expect_equal(as.numeric(lh["linf", "M"]), 78)
  r <- ss3_resid(ss3ank)
  for(ty in c("bubble", "len", "qq"))
    expect_silent(suppressWarnings(ggplot2::ggplot_build(ss3plot_resid(r, type = ty, lhpar = lh))))
  skip_no_data(); e <- load_ss3("pac.hke")
  expect_silent(ggplot2::ggplot_build(ss3plot_resid(ss3_resid(e$ss3phk, quant = "age"))))
})
