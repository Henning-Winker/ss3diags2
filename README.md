# ss3diags2

**Stock Synthesis diagnostics and data extraction in FLR**

`ss3diags2` is the next generation of [`ss3diags`](https://github.com/JABBAmodel/ss3diags). It extracts observations
and model outputs from Stock Synthesis (SS3) reports into long `data.table` tables that convert directly to
FLR classes, shaped to match the inputs of [`JABBA3`](https://github.com/Henning-Winker/JABBA3) and
[`FLicc`](https://github.com/Henning-Winker/FLicc). Diagnostics in `ggplot2` follow.

## Installation

```r
remotes::install_github("Henning-Winker/ss3diags2")
```

## Example data

- `ss3simple`: the simulated two-sex "simple_small" model of r4ss (SS 3.30.26; one fishery, two surveys, length and age
  compositions). It is the example of the PIFSC ss3diags, so results can be compared with ss3diags; the user guide
  runs on it.
- Northern shrimp pra.27.3a4a (SS 3.30.22): four seasons, protandric hermaphrodite (two sexes), two areas, six
  area-specific fleets with retention. Report and SS3 input files:
  `readRDS(system.file("extdata", "pra.27.3a4a", "pra.27.3a4a.rds", package = "ss3diags2"))`; user guide section
  "Season-, sex- and area-structured models" (F-at-age, partial F, F by area, stock by area, partial impact).
- Cookbook examples (Carvalho et al. 2021, copied from ss3diags): North Atlantic shortfin mako (`ss3sma`,
  `retro.sma`, `aspm.sma`) and Pacific hake (`ss3phk`, `retro.phk`, `aspm.phk`), to replicate the paper.

## Observations from an SS3 report

```r
library(ss3diags2)
data(ss3sma)                          # r4ss::SS_output() of ICCAT N Atlantic shortfin mako (cookbook)

ss3_fleets(ss3sma)                      # fleet table: type, area, units, data available

catch <- ss3_catch(ss3sma)              # FLQuants by fleet
idx   <- ss3_indices(ss3sma)            # FLIndices of FLIndexBiomass (index, index.var, index.q, timing)
lfd   <- ss3_lencomp(ss3sma)            # FLQuants by fleet: len x year x unit (sex) x season x area

data(ss3phk); data(retro.phk)
afd   <- ss3_agecomp(ss3phk, scale = "prop")

# every extractor can return the long table behind the FLR object
head(ss3_lencomp(ss3sma, long = TRUE))
```

## Ready-made inputs for JABBA3 and FLicc

```r
inp <- ss2jabba3(ss3sma)                # annual catch by fleet (weight) + FLIndices
# fit <- JABBA3::jb3_fit(JABBA3::jb3_data(inp$catch, inp$indices, stock = "sma"), ...)

lb <- ss2flicc(ss3sma)                  # lfd by gear (common bins and years) + catch_by_gear
# stk <- FLicc::stocklen(lb$lfd, lhpar)
```

## Model outputs as FLStockR

```r
stk  <- ss3stk(ss3sma)                 # MLE with log-scale SEs; refpts Ftgt, Btgt, MSY, B0, R0
stki <- ss3stk(ss3sma, output = "iters")  # joint MVLN draws (ss3mvln, as ss3om::ssmvln)
```

## Plots (ggplot2, FLR inputs, r4ss colours)

```r
ss3plot_data(catch = catch, indices = idx, lencomp = lfd)
ss3plot_catch(catch, by = "season")                      # by any of "unit", "season", "area"
ss3plot_indices(idx, fit = ss3_indices(ss3sma, what = "exp"))
ss3plot_comps(lfd, by = "unit", type = "pooled")         # "annual", "bubble", "mean"
ss3plot_trj(stki)                                        # also list(SS3 = stk, JABBA3 = jb3stk)
ss3plot_kobe(stki)
ss3plot_kobe(ss3stk(ss3sma, output = "iters", fwd = TRUE), lag = 1)   # F(y) vs B(y + 1)
ss3plot_production(stk, type = "phase")                  # or "dyn"; Pella-Tomlinson fit overlaid
```

## Production function and JABBA3 priors

```r
data(ss3phk); data(retro.phk)
sp <- ss3_spm(ss3phk)        # Pella-Tomlinson fit to the SS3 equilibrium curve, checked against SS3 MSY, Bmsy, B0
sp$pars                      # r, K, m, MSY, Bmsy, BmsyK, Hmsy
sp$priors                    # c(median, log.sd, on) for JABBA3::jb3_priors(r = , K = , m = )
```

See [DESIGN.md](DESIGN.md) for the data model, decisions and roadmap (life-history `FLPar`, MLE `FLStockR`, MVLN
without dplyr, ggplot2 diagnostics, OSA residuals and posterior predictive checks).

## User guide

The user guide (<https://henning-winker.github.io/ss3diags2/>; source `inst/guide/ss3diags2_user_guide.Rmd`, R code with
comments in `inst/ss3diags2_user_guide.R`, both updated by `inst/guide/updguide.R`) follows the JABBA3
guide: Part I goes through the diagnostics with checks (data, fits, OSA residuals, runs tests, retrospective,
hindcast; cookbook examples), Part II covers stock status and the reference point basis, the timing of F and SSB in
Kobe plots (`ss3plot_kobe(lag = 1)`), the production function, season-, sex- and area-structured models (F by area,
stock by area, partial impact of fleets after Wang et al. 2009) and the FLR outputs.

## Acknowledgements

The joint MVLN (`ss3mvln()`) follows Kell et al. (2023, Collect. Vol. Sci. Pap. ICCAT 80(4): 837-858) and
`ss3om::ssmvln()`; the year-by-year option is the delta-MVLN of Walter and Winker (2020, ICCAT SCRS/2019/145) and
Winker et al. (2019, ICCAT SCRS/2019/093), implemented in `ss3diags::SSdeltaMVLN()`.

## License

EUPL (>= 1.2)
