# ss3diags2 design notes (draft, 11 Nov 2026)

Working notes for the structure of ss3diags2, the FLR/ggplot2 successor of
[ss3diags](https://github.com/JABBAmodel/ss3diags). Items marked **[open]** need a decision before the code that
depends on them is written.

## Goals
1. Extract observations, then model outputs, from Stock Synthesis reports (`r4ss::SS_output()`, here `ss3rep`) into
   FLR objects, so an SS3 model can be refitted with JABBA3 or FLicc and compared on the same data.
2. One data model: long `data.table` tables whose columns are the FLQuant dimensions
   (`len`/`age`, `year`, `unit`, `season`, `area`, `iter`, `data`) plus SS3 bookkeeping columns (`fleet`, `name`,
   `obs`, `exp`, `se`, `nsamp`, ...). Conversion to `FLQuant(s)`, `FLIndices` and `FLStock` is then a reshape.
3. Diagnostics (runs test, joint residuals, retro, hindcast/MASE, MVLN, Kobe) in ggplot2, harmonised in names and
   style with JABBA3 and FLicc; OSA residuals and posterior predictive checks added.
4. No dplyr. Base R in the core, data.table for the large tables (comps, MVLN/MCMC draws).

## Conventions (as JABBA3 and FLicc)
- FLR S4 classes in and out; plain lists with an S3 tag for fit-like results.
- Function prefixes: `ss3_` extract from an SS3 report; `ss2<model>()` shape the extraction as inputs of another
  model; `ssplot_` plots (to come). Every extractor has `long = TRUE` to return the underlying data.table.
- FLR dimensions: sex -> `unit` ("F", "M"; "unique" for combined-sex data or single-sex models), SS3 season ->
  `season` ("all" when one season), SS3 area -> `area` ("unique" when one area), draws -> `iter`, fleets -> elements
  of `FLQuants` / `FLIndices` (as FLicc gears and JABBA3 catch fleets).
- Units: catch "t" (biomass) or "1000" (numbers, SS3 convention); comps "numbers" or "prop"; index units from SS3
  survey units ("numbers", "biomass", ...); index SEs "se" (log scale).

## Observations (done, v0.1)

| ss3diags2 | SS3 source (r4ss) | Returns | JABBA3 input | FLicc input |
|---|---|---|---|---|
| `ss3_fleets()` | `FleetNames`, `fleet_type`, `catch_units`, `survey_units`, data tables | data.table: fleet, name, type, area, units, n obs by data type | choose catch / index fleets | choose gears (catch fleets with lengths) |
| `ss3_catch()` | `catch` (Obs, Exp, ret/kill/sel bio and num, F) | FLQuants by fleet (quant "all", season, area) | `catch` (FLQuants by fleet) | `catch_by_gear` (FLQuants) |
| `ss3_indices()` | `cpue` (Obs, Exp, SE, SE_input, Calc_Q, Month, Use) | FLIndices of FLIndexBiomass: index, index.var, index.q, range startf/endf | `indices` (FLIndexBiomass, numbers or weight, SE in `index.var`) | (none) |
| `ss3_lencomp()` | `lendbase` | FLQuants by fleet: len x year x unit(sex) x season x area | (Selex, to come) | `lfd` (FLQuants by gear, len x year) |
| `ss3_agecomp()` | `agedbase` (marginal ages) | FLQuants by fleet: age x year x unit x season x area | (none) | (none) |
| `ss2jabba3()` | the above | list(catch, indices, catch.total, fleets), annual, catch in weight | ready for `jb3_data(catch, indices)` | (none) |
| `ss2flicc()` | the above | list(lfd, catch_by_gear, years): sexes, seasons, partitions summed; common bins and years | (none) | ready for `stocklen(lfd, lhpar)` and fitting |

Decisions (11 Nov):
- Index SE: `index.var` holds the **input** SE (`SE_input`) by default, as JABBA3 estimates an additive `sdi` on top
  (`se = "total"` gives the SE used in the SS3 likelihood). Older r4ss outputs lack `SE_input`; then the total SE is
  used, with a message.
- Index timing: `range[c("startf", "endf")]` = (month - 1)/12 (SS3 3.30 month), else time - year.
- Observations excluded in SS3 (`Use` < 0) are dropped from FLIndices by default (`excluded = "keep"` keeps them).
- Compositions: SS3 proportions x sample size (`nsamp = "adj"`, after variance adjustment, or `"input"`) give
  numbers (`scale = "n"`, default; FLicc fits counts). Combining sexes or several samples in one cell sums numbers
  and recomputes proportions from them.
- Partitions (discard, retained) give separate elements `fleet.discard`, `fleet.retained`. Ageing-error definitions
  stay in the long table (`ageerr`): they are often year-specific (Pacific hake has one per year).
- Catch for JABBA3 and FLicc must be in weight: with `catch = "auto"` fleets reporting numbers use the SS3 retained
  catch in weight (`ret_bio`), with a message.
- FLicc fits every gear in every year: `ss2flicc(years = "common")` keeps years in which all gears have a sample;
  `"all"` keeps every year with `NA` for gears without a sample **[open: check FLicc treats NA years as missing]**.

## Model outputs (done, v0.2)

| ss3diags2 | SS3 source | Returns | Matches |
|---|---|---|---|
| `ss3mvln()` | `derived_quants`, `CoVar`, `catch` | list(kb, mle: data.tables; refpts: FLPar; basis) | `ss3diags::SSdeltaMVLN()` (no dplyr) |
| `ss3stk()` | the above | FLStockR (age "1"): ssb = SSB, fbar = F, rec, catch; refpts Ftgt, Btgt, MSY, B0, R0 (+ Blim, Bpa) | `FLRef::ss2FLStockR()`, `JABBA3::jb3_FLStockR()` |
| `ss3_eqcurve()` | `equil_yield` | data.table F, SSB, yield (+ Fmsy, Bmsy, MSY, B0) | production function |

- MLE: `iters = 0`, log-scale SEs of SSB, F, recruitment, B/Btgt and F/Ftgt in `attr(x, "se.log")` (Wald intervals).
- MVLN: `iters = n`, draws as iterations, MLE in `attr(x, "mle")`, jointly sampled B/Btgt and F/Ftgt in
  `attr(x, "status")` (SSB and F draws are marginal, so SSB/Btgt from them would lose the B-F correlation).
- MLE values and reference points equal `SSdeltaMVLN()` to machine precision (3 test models). The draws differ slightly:
  `SSdeltaMVLN()` subsets CoVar with `paste0(c("Bratio", "F"), "_", yrs)`, which recycles element-wise (Bratio in odd,
  F in even positions), so its neighbouring-year correlations use a subset; `ss3mvln()` uses all years
  (terminal-year medians agree within 0.1%).

## Plots (done, v0.2): `ss3plot_` on FLR objects

All plots take FLR objects (not the SS3 report), so they also work for JABBA3, FLicc and FLRef objects. `by` keeps
FLQuant dimensions `"unit"`, `"season"`, `"area"` separate (other dimensions are summed). Colours follow r4ss:
fleets `rich.colors.short` (`ss3cols()`), seasons black / blue4, green3, orange2, red3 (`ss3seascols()`), unit
unique green, F red, M blue (`ss3unitcols`), observed compositions grey bars, fitted index blue. Theme `theme_ss3()`
= JABBA3 theme. Stock plots mirror JABBA3 (`jbplot_trj`, `jbplot_kobe`, `jbplot_spphase`, `jbplot_spdyn`).

| Plot | Input | Options |
|---|---|---|
| `ss3plot_data()` | catch FLQuants, FLIndices, len/age FLQuants | `by` (panels) |
| `ss3plot_catch()` | FLQuants by fleet | `type` bar / line / share, `by` (panel rows) |
| `ss3plot_indices()` | FLIndices (+ fit FLIndices) | CIs from index.var, `by` (colour), `combine`, `log`, `time` |
| `ss3plot_comps()` | FLQuants of numbers (+ fit) | `type` pooled / annual / bubble / mean, `by` (colour) |
| `ss3plot_nsamp()` | FLQuants of numbers | `by` |
| `ss3plot_trj()` | FLStockR, FLStocks, list | SSB, F, Rec, Catch, B/Btgt, F/Ftgt, B/B0, Catch/MSY; nested 90/60% bands; runs compared |
| `ss3plot_kobe()` | FLStockR with iters | points or filled HDR contours, quadrant percentages |
| `ss3plot_production()` | FLStockR (+ equilibrium curve) | `type` phase / dyn, `rel` |

## Production function (done, v0.3)

- `ss3_eqcurve()` picks one SS3 profile: `equil_yield` stacks repeated F grids (SPR loops 0-2), single benchmark rows
  (3-6) and an extension from Fmsy to high F (7, some versions). Default: the last grid from F = 0 (converged; it
  reproduces the SS3 MSY) plus extension loops. Loop 0 can differ (mako: another F scale, yield 33 at Fmsy).
- `ss3_spm()` fits Pella-Tomlinson (r, K, m; Fox at m = 1) to yield against SSB, as `FLRef::asem2spm()`; extrapolates to
  SSB = 0 and K; checks MSY, Bmsy, B0, Bmsy/B0 against SS3 (warning above 5%); returns JABBA3 priors
  `c(median, log.sd, on)` for r, K (CVs from SS3 SEs of MSY and B0) and m (off by default).
- Results: hake R2 0.99993, benchmarks within 0.2%; mako R2 0.9999, within 0.6% (SSB in numbers, so r and K priors are
  off, the shape carries over); herring R2 0.9997 but the equilibrium table implies unfished SSB 1.62 x SSB_unfished
  and MSY 13% above Dead_Catch_MSY: the table and the SS3 MSY calculation use different settings **[open: which
  forecast.ss benchmark years?]**.
- Caveat: with m < 1 the curve has an infinite slope at the origin, so extrapolation to very low SSB is optimistic.
- `ss3plot_production(spm = TRUE)` overlays the fit (dashed) on the SS3 curve.
- **[open]** production function in exploitable biomass (rebuild equilibrium numbers at age per F from `endgrowth`,
  `ageselex`, M and F shares; or FLBRP + `asem2spm(quant = "vb")`), so yield and biomass share units (mako) and JABBA3
  gets EB/SB mapping for Selex; realised SP(t) fit with process error as a diagnostic.

## Example data: ss3ank (v0.4)

Black anglerfish ank.27.78abd, SS3 3.30.22, 1950-2025. Checked against the input files (data.ss: Nseas = 1 with 4
sub-seasons, Ngenders = 2, Nareas = 1):
- not seasonal; sex structure only in the survey length compositions (Sexes = 3); fishery compositions are
  sex-combined with discard and retained partitions; observed discards (`ss3_discards()`); input SEs present.
- the report has no CoVar, derived quantities are 0 and `equil_yield` is empty (run without hessian/sdreport):
  `ss3stk()` falls back to the time series (SSB, recruitment, dead catch, F = sum of apical F, B0, R0);
  no MVLN, Kobe or production function. **[open: rerun with hessian; a seasonal variant would test the season code]**
- triggered fixes: r4ss >= 1.50 names (`dead_bio`), compositions completed on all data bins (SS3 tail compression
  drops empty bins), fit lines broken at gaps.

## Composition residuals (v0.5)

- `ss3_resid()`: OSA (default) or Pearson residuals of length / age compositions as FLQuants by fleet and partition.
  OSA as FLicc (Trijoulet et al. 2023): conditional binomial (multinomial) or beta-binomial (Dirichlet-multinomial)
  along increasing bins, females before males in joint-sex compositions; DM precision from the SS3 effective N,
  a0 = N (Neff - 1)/(N - Neff) (`Nsamp_DM`, `DM_effN`), equal to theta x N for the linear parameterisation (checked on
  ss3ank). Counts round(obs x N). Validated on simulated compositions (mean 0, sd 1).
- `ss3plot_resid()`: bubbles (panels per fleet x unit with data, years with data), cohort tracks from `ss3_lhpar()`
  by unit (F solid, M dashed; diagonals for ages), residuals pooled over years by bin (Stouffer), QQ.
- `ss3_lhpar()`: FLPar linf, k, t0, L50, mat.slope, a, b, M, CVmin, CVmax, A1, A2, L1, L2 (unit dimension with two
  sexes); VB approximated from `endgrowth` for other growth models; missing male parameters set to female values.
- ss3ank: OSA means 0.1-0.4 and survey residuals negative at 10-25 cm, positive above 30 cm in both sexes: fewer small
  fish than expected (conditional reading), a selectivity or recruitment signal; hake and mako residuals ~N(0,1).
- **[open]** OSA for conditional age-at-length; OSA for indices (needs the SS3 likelihood: lognormal with total SE).

## MVLN and FLStockR (v0.6)

- `ss3mvln()` (authors L.T. Kell, H. Winker), data.table, no dplyr:
  - `method = "joint"` (default), the `ss3om::ssmvln()` approach: one MVLN over SSB, recruitment and F of all years and
    the reference points (Btgt, Ftgt, MSY, B0, R0) with the full CoVar correlation matrix (CoVar holds one
    triangle; filled symmetric; nearest positive definite), log-scale SD sqrt(log(1 + CV^2)), median at the MLE.
    Iterations are coherent trajectories, reference points have iterations, B/Btgt and F/Ftgt are ratios within
    iterations. Extends ss3om (Recr, Bratio, F only) with SSB and the reference points; matches it on the shared
    quantities (sma 2015: Bratio median 0.961 vs 0.965, CV 0.0486 vs 0.0487; F/Fmsy 4.28 vs 4.25).
  - `method = "annual"`: `SSdeltaMVLN()` port (year-by-year).
  - Finding: SSdeltaMVLN attenuates the B-F correlation (sma 2015: -0.24 vs CoVar -0.86; hake -0.26 vs -0.76),
    the joint method keeps it (-0.86). Likely causes: the covariance entries [2,3] and [1,3] are swapped (cov(B, Fref)
    and cov(F, Fref)), and the recycled `paste0(status, "_", yrs)` subset of correlations. **[open: fix SSdeltaMVLN?]**
- `ss3stk()` (as `FLRef::ss2FLStockR()`): FLStockR from an SS3 report (`output = "mle"` from derived quantities with
  SEs, time-series fallback; `"iters"` from `ss3mvln()`) or from an `ss3mvln` object; `rel`, `thin`, Blim/Bpa.
  Replaces `ss3_FLStockR()`.
- `ss3_flqs()`: SSB, F, Rec, Catch, B/Btgt, F/Ftgt, B/B0, Catch/MSY as FLQuants (iterations kept, ratios within
  iterations) for ggplotFL `plot()` and `ss3plot_trj()`; used by the stock plots.

- `fwd = TRUE` (ss3mvln, ss3stk) adds the SS3 forecast years; `attr(stk, "endyr")` marks the last assessment year.
- `ss3plot_kobe(lag = 1)`: F(y) against B/Btgt(y + 1), the start-of-year biomass after the fishing in year y; the cloud
  pairs F(Y) and B(Y + 1) within iterations (joint MVLN; warning for SSdeltaMVLN draws). sma 2015: P(red) 78% at
  lag 0, 93% at lag 1.

## Retro, hindcast, runs tests (v0.7)

- `ss3_retro()`: FLStocks of MLE `ss3stk()` per peel (end year from `Retro_year`; `fwd = TRUE` for the year after).
- `ss3_mohn()`: Mohn's rho and forecast rho for any `ss3_flqs()` quantity; equals `SShcbias()` (sma, hake). Note:
  SShcbias "F" is the reported F (F/Fmsy for sma), `ss3_mohn("F")` is absolute F; use "F/Ftgt" for the ratio.
- `ss3_hindcast()` + `ss3plot_hcxval()`: one-year-ahead predictions of indices or mean length/age, MASE and
  MASE.adj; equals `SSmase()`.
- `ss3_iresid()`, `ss3_runs()`, `ss3plot_runs()`, `ss3plot_jres()`: log residuals of indices or mean length/age,
  runs test (own Wald-Wolfowitz, no randtests) with 3-sigma limits, joint residuals with RMSE; equal `SSrunstest()`.

## Season-, sex- and area-structured models: ss3pra (v0.8, 12 Nov)

Test case: northern shrimp pra.27.3a4a (`data(ss3pra)`, SS 3.30.22): 4 seasons (spawning in season 1), 2 sexes
(protandric hermaphrodite, born male), 2 areas with time-varying recruitment distribution, 6 area-specific catch
fleets (DK/NOR/SWE x 3a/4a) with retention, 2 surveys with sex-specific length compositions. F basis: Fbar 1-3 as
Z - M; MSY proxy F at B30%.

Structure map (what each dimension means and where it goes):

| SS3                    | FLR dimension | extraction                              | note |
|------------------------|---------------|-----------------------------------------|------|
| Sexes 0 / 1-2 (data)   | unit          | unique / F, M                           | combined-sex data stay "unique" |
| Seas                   | season        | catch, comps, fa, sel (age)             | seasons are not years (no step2year needed) |
| Area                   | area          | catch, fa, byarea; fleet area from `fleet_area` | fleets are area-specific, so observations need no area split |
| Part 0/1/2             | fleet name    | `fleet.all`, `fleet.discard`, `fleet.retained` | |
| Morph / platoon        | (averaged)    | fa, sel                                 | one morph per sex here |
| Era INIT/TIME/FORE     | year range    | `era =` in `ss3_fa()`                   | |

All existing extractors and plots ran on ss3pra without changes (catch by season x area, indices, length comps by
unit/season/part, OSA residuals with growth tracks for both sexes, ss3stk/ss3mvln, Kobe with lag, production
function, runs tests). The SPM fit warns (Bmsy > 5% from the SS3 benchmark; MSY proxy is F at B30%).

New (from the r4ss HTML tabs and the GSA 17-18 hake advice):
- `ss3_fa()`: F-at-age by fleet as FLQuants (age, year, unit, season, area), from `fatage` (getFa of the hake Rmd);
  `annual = TRUE` sums seasons weighted by duration. `ss3plot_fa()`: F-at-age by fleet and season.
- `ss3_partialF()` + `ss3plot_partialF()`: partial F by fleet, gear (named list) or area. `method = "catch"` splits the
  reported F (exploitation$annual_F) by catch numbers over the F ages, so partial F's add up to F (tested);
  `"fbar"` / `"apical"` give the F experienced in the area of the fleet (F by area with `groups = "area"`), optionally
  by season.
- `ss3_impact()` + `ss3plot_impact()`: partial impact after Wang et al. (2009) without re-running SS3 or FLasher.
  Annual Z from the SS3 numbers at age (seasons, areas, sexes pooled), SS3 M at age (season-weighted), fishing
  Z - M split among fleets by catch at age in numbers, SS3 spawning output at age at spawning time, BH/Ricker SR with
  the SS3 deviations kept. Validation on ss3pra: no fishing (init = "virgin") reproduces `Dynamic_Bzero` to 1e-5;
  all fleets reproduce the SS3 SSB (exactly with init = "model"). Approximation: movement, sex change and in-year
  timing held at their fitted values. `type = "remove"` for one-group-out.
- `ss3_sel()` + `ss3plot_sel()`: selectivity, retention, discard mortality at length (time blocks) and selectivity at
  age, by unit. `ss3_byarea()` + `ss3plot_byarea()`: SSB, biomass, recruitment by area (FLQuants with area dimension).
  `ss3plot_sr()`: stock-recruit curve and recruitment deviations.
- `ss3plot_catch(by = c("season", "area"))`: two dimensions now give a grid.

Still missing from the r4ss tabs: numbers at age / mean age by area (Numbers), hermaphroditism transition and
maturity/fecundity at age (Bio), SPR curves (SPR), parameter table with priors (Pars), discard fits plot (Discard).

## Flagged for later: NUTS (11 Nov)

- `ss3nuts(dir, ...)`: wrapper of `adnuts::sample_nuts()` on the SS3 executable (after an `-hbf` MLE run), chains,
  warmup, metric = "mle"; adnuts and rstan in Suggests only (`.need()`), no hard dependency.
- `ss3_posteriors(dir)` / from an adnuts fit: read `derived_posteriors.sso` (as `ss3diags::SSreadMCMC()`) into the
  `ss3mvln` draws table (method = "nuts": SSB, Recr, F, Bratio by year and the reference points per iteration), so
  `ss3stk()`, `ss3plot_trj()` and `ss3plot_kobe()` (incl. `lag = 1`) work unchanged and NUTS can be compared with
  the joint MVLN. Test data: `ss3diags::gob.her` (`mcmc.her`, derived posteriors).
- NUTS diagnostics (trace, rank, R-hat, ESS, divergences) in the style of `JABBA3::jbplot_nuts()`.
- Blocker in this workspace: no SS3/ADMB binary (release assets not reachable); test where SS3 runs.

## Next (in this order)
1. **Life history `FLPar`** (`ss3_lhpar()`), from `Growth_Parameters`, `MGparmAdj` and `endgrowth`: linf, k, t0,
   L50 and slope (or a50), M (constant or by age), a, b, CVL (CVmin/CVmax), A1/A2/L1/L2, by sex in a `unit`
   dimension. For growth other than von Bertalanffy (Richards, age-specific K, empirical length at age), a helper
   that approximates VB by fitting mean length at age from `endgrowth` and reports the fit **[open: weighting by
   numbers at age?]**. Target: FLicc `lhpar` and JABBA3 life-history priors.
2. **Conditional age-at-length** (`condbase`) and size compositions; mean length/age (Francis TA1.8) as FLQuants.
3. ~~Model outputs, MLE, FLStockR~~ (done).
4. **MVLN** (done: `ss3mvln()` joint and annual): remaining: from `ss3om::ssmvln()` / `mvlnFLQuants()` (data.table) and
   `FLRef::ssmvln()` (the SSdeltaMVLN port), returning draws as FLQuants with iters and the FLStockR with iters.
   Must work for seasonal ("interannual", seasons-as-years) models, i.e. map SS3 time steps to year/season as
   `ss3diags::SSdiagsStep2Year()` does.
5. **Age-structured FLStock** where the model allows (from `ss3om::buildFLSss330()` logic: unit = sex x birth
   season, season, area).
6. **Diagnostics in ggplot2** on the long tables: runs test, joint residuals, retro (Mohn's rho), hindcast (MASE),
   model comparison, ensemble, Kobe; then OSA residuals and posterior predictive checks (full posteriors only, as
   JABBA3). Harmonise with `jbplot_*` (JABBA3) and `plot_*` (FLicc) **[open: prefix `ssplot_` or `plot_ss3_`?]**.

## Sources to reuse
- `ss3om` (flr): data.table readers, `getDimnames()`, `codeUnit()` (unit from sex and morph), `ssmvln()`,
  `mvlnFLQuants()`.
- `FLRef`: `ss2FLStockR()`, `ssmvln()`, `ss2ices()`.
- `ss3diags` (JABBAmodel): diagnostic statistics (`ssruns_sig3()`, `SSmase()`, `SShcbias()`, `SScompsTA1.8()`).
