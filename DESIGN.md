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
| `ss3_mvln()` | `derived_quants`, `CoVar`, `catch` | list(kb, mle: data.tables; refpts: FLPar; basis) | `ss3diags::SSdeltaMVLN()` (no dplyr) |
| `ss3_FLStockR()` | the above | FLStockR (age "1"): ssb = SSB, fbar = F, rec, catch; refpts Ftgt, Btgt, MSY, B0, R0 (+ Blim, Bpa) | `FLRef::ss2FLStockR()`, `JABBA3::jb3_FLStockR()` |
| `ss3_eqcurve()` | `equil_yield` | data.table F, SSB, yield (+ Fmsy, Bmsy, MSY, B0) | production function |

- MLE: `iters = 0`, log-scale SEs of SSB, F, recruitment, B/Btgt and F/Ftgt in `attr(x, "se.log")` (Wald intervals).
- MVLN: `iters = n`, draws as iterations, MLE in `attr(x, "mle")`, jointly sampled B/Btgt and F/Ftgt in
  `attr(x, "status")` (SSB and F draws are marginal, so SSB/Btgt from them would lose the B-F correlation).
- MLE values and reference points equal `SSdeltaMVLN()` to machine precision (3 test models). The draws differ slightly:
  `SSdeltaMVLN()` subsets CoVar with `paste0(c("Bratio", "F"), "_", yrs)`, which recycles element-wise (Bratio in odd,
  F in even positions), so its neighbouring-year correlations use a subset; `ss3_mvln()` uses all years
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

## Next (in this order)
1. **Life history `FLPar`** (`ss3_lhpar()`), from `Growth_Parameters`, `MGparmAdj` and `endgrowth`: linf, k, t0,
   L50 and slope (or a50), M (constant or by age), a, b, CVL (CVmin/CVmax), A1/A2/L1/L2, by sex in a `unit`
   dimension. For growth other than von Bertalanffy (Richards, age-specific K, empirical length at age), a helper
   that approximates VB by fitting mean length at age from `endgrowth` and reports the fit **[open: weighting by
   numbers at age?]**. Target: FLicc `lhpar` and JABBA3 life-history priors.
2. **Conditional age-at-length** (`condbase`) and size compositions; mean length/age (Francis TA1.8) as FLQuants.
3. ~~Model outputs, MLE, FLStockR~~ (done).
4. **MVLN** (done for annual models): remaining: `ssmvln()` from `ss3om::ssmvln()` / `mvlnFLQuants()` (data.table) and
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
