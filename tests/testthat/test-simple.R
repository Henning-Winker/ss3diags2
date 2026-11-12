# simulated two-sex model (r4ss simple_small), the PIFSC ss3diags example
data(ss3simple, package = "ss3diags2")

test_that("ss3simple runs through observations, residuals and status", {
  expect_equal(ss3_fleets(ss3simple)$name, c("FISHERY", "SURVEY1", "SURVEY2"))
  expect_length(ss3_indices(ss3simple), 2)
  expect_true(all(c("F", "M") %in% dimnames(ss3_agecomp(ss3simple)[[1]])$unit))
  expect_s4_class(ss3_resid(ss3simple, "age"), "FLQuants")
  expect_equal(nrow(ss3_runs(ss3simple)), 2)
  stk <- ss3stk(ss3simple)
  expect_s4_class(stk, "FLStockR")
})

test_that("impact with initial F: all fleets = SS3 SSB, no fishing = SS3 dynamic B0", {
  imp <- ss3_impact(ss3simple)
  expect_equal(imp$check$projected, imp$check$ss3, tolerance = 1e-8)
  expect_equal(imp$check$dynB0, imp$check$ss3_dynB0, tolerance = 1e-4)
  pf <- ss3_partialF(ss3simple)
  expect_equal(c(pf[[1]]), c(attr(pf, "total")), tolerance = 1e-8)
})
