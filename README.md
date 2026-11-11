# ss3diags2

**Stock Synthesis diagnostics and data extraction in FLR**

`ss3diags2` is the next generation of [`ss3diags`](https://github.com/JABBAmodel/ss3diags). It extracts observations
and (soon) model outputs from Stock Synthesis (SS3) reports into long `data.table` tables that convert directly to
FLR classes, shaped to match the inputs of [`JABBA3`](https://github.com/Henning-Winker/JABBA3) and
[`FLicc`](https://github.com/Henning-Winker/FLicc). Diagnostics in `ggplot2` follow.

## Installation

```r
remotes::install_github("Henning-Winker/ss3diags2")
```

## Observations from an SS3 report

```r
library(ss3diags2)
data(natl.sma, package = "ss3diags")   # ss3sma: r4ss::SS_output() of ICCAT N Atlantic shortfin mako

ss3_fleets(ss3sma)                      # fleet table: type, area, units, data available

catch <- ss3_catch(ss3sma)              # FLQuants by fleet
idx   <- ss3_indices(ss3sma)            # FLIndices of FLIndexBiomass (index, index.var, index.q, timing)
lfd   <- ss3_lencomp(ss3sma)            # FLQuants by fleet: len x year x unit (sex) x season x area

data(pac.hke, package = "ss3diags")
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

See [DESIGN.md](DESIGN.md) for the data model, decisions and roadmap (life-history `FLPar`, MLE `FLStockR`, MVLN
without dplyr, ggplot2 diagnostics, OSA residuals and posterior predictive checks).

## License

EUPL (>= 1.2)
