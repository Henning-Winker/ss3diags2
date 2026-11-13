# synthetic conditional age-at-length data on ss3simple: age given length from von Bertalanffy growth with
# normal length-at-age and exponential survival; observations drawn from the model (multinomial, N per bin)
fake_caal <- function(r, fleet = 2, years = 2013:2022, N = 25, seed = 3){
  set.seed(seed)
  lh <- ss3_lhpar(r); a <- lh@.Data
  linf <- a["linf", 1, 1]; k <- a["k", 1, 1]; t0 <- a["t0", 1, 1]; M <- a["M", 1, 1]
  ages <- 1:15; lb <- seq(20, 80, by = 5)
  La <- linf*(1 - exp(-k*(ages + 0.5 - t0)))
  db <- data.table::rbindlist(lapply(years, function(y) data.table::rbindlist(lapply(lb, function(l){
    p <- exp(-M*ages)*stats::dnorm(l + 2.5, La, 0.1*La); p <- p/sum(p)
    o <- as.vector(stats::rmultinom(1, N, p))/N
    data.table::data.table(Yr = y, Month = 7, Seas = 1, Fleet = fleet, Area = 1, Repl. = 1, Sexes = 0, Sex = 1,
                           Part = 0, Ageerr = 1, Lbin_lo = l, Lbin_hi = l + 5, Bin = ages, Obs = o, Exp = p,
                           Pearson = (o - p)/sqrt(p*(1 - p)/N), Nsamp_adj = N, Like = 1)
  }))))
  r$condbase <- as.data.frame(db)
  r
}
