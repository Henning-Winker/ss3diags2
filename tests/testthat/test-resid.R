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
  data(ss3phk, package = "ss3diags2", envir = environment())
  d <- ss3_resid(ss3phk, quant = "age", long = TRUE)
  ab <- ss3phk$agedbase; s <- ab[ab$Fleet == 1 & ab$Yr == 2010, ][1, ]
  th <- exp(ss3phk$parameters["ln(EffN_mult)_1", "Value"])
  expect_equal(d[name == "Fishery" & year == 2010]$a0[1], th*s$Nsamp_adj, tolerance = 0.01)
  r <- ss3_resid(ss3phk, quant = "age")
  expect_s4_class(r, "FLQuants")
  expect_equal(attr(r, "type"), "osa")
  pe <- ss3_resid(ss3phk, quant = "age", type = "pearson", long = TRUE)
  expect_equal(sum(is.finite(pe$data)), sum(is.finite(ab$Pearson)))
})

test_that("ss3_lhpar returns von Bertalanffy parameters by unit and residual plots build", {
  data(ss3simple, package = "ss3diags2", envir = environment())
  lh <- ss3_lhpar(ss3simple)
  expect_equal(dimnames(lh)$unit, c("F", "M"))
  r <- ss3_resid(ss3simple)
  expect_equal(dimnames(r[["SURVEY1"]])$unit, c("F", "M"))
  for(ty in c("bubble", "len", "qq"))
    expect_silent(suppressWarnings(ggplot2::ggplot_build(ss3plot_resid(r, type = ty, lhpar = lh))))
})
