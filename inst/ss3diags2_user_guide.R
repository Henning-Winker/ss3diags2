# ss3diags2 user guide: R code with comments
# Generated from inst/guide/chapters by inst/guide/updguide.R; edit the chapters, not this file.
# Rendered guide: https://henning-winker.github.io/ss3diags2/

library(ss3diags2)
library(ggplot2)

####################################################################################################
# SS3DIAGS2: STOCK SYNTHESIS DIAGNOSTICS IN FLR AND GGPLOT2
####################################################################################################

# ss3diags2 is the successor of ss3diags (Winker et al.; Carvalho et al. 2021). It reads the
# output of Stock Synthesis (SS3; Methot and Wetzel 2013) as returned by r4ss::SS_output() and
# provides:

# * Observations as FLR objects. Catch, discards and length and age compositions as FLQuants by
#   fleet, indices as FLIndices, life-history parameters as FLPar, all with the SS3 structure in
#   the FLQuant dimensions (sex as unit, season, area). The same objects are the inputs of JABBA3
#   and FLicc.

# * Residual diagnostics. One-step-ahead (OSA) residuals for compositions (Trijoulet et al.
#   2023), runs tests, joint residual plots, retrospective analysis (Mohn's rho) and hindcast
#   cross-validation (MASE), as in the diagnostics cookbook (Carvalho et al. 2021).

# * Stock status with uncertainty. FLStockR objects with reference points, maximum likelihood
#   estimates or draws from a multivariate lognormal (MVLN) approximation of the SS3 Hessian (Kell
#   et al. 2023; Walter and Winker 2020), trajectories, Kobe plots and the production function.

# * Fishing mortality and impacts. F-at-age by fleet, partial F by fleet, gear or area, F and
#   stock by area, and the partial impact of fleets on the spawning stock (Wang et al. 2009), from
#   the SS3 report without re-running SS3.

# * data.table and ggplot2. Long tables with FLQuant dimension columns (len/age, year, unit,
#   season, area, iter, data) that convert directly to FLR; no dplyr.

# How ss3diags functions map to ss3diags2:

# | ss3diags | ss3diags2 |
# | SSplotRunstest(), SSrunstest() | ss3_runs(), ss3plot_runs() (indices, mean length, mean age) |
# | SSplotJABBAres(), SSrmse() | ss3plot_jres() |
# | SSplotRetro(), SShcbias() | ss3_retro(), ss3_mohn(), ss3plot_retro() |
# | SSplotHCxval(), SSmase() | ss3_hindcast(), ss3plot_hcxval() |
# | SSdeltaMVLN(), SSplotKobe() | ss3mvln(), ss3stk(), ss3plot_trj(), ss3plot_kobe(), ss3plot() |
# | SSsettingsBratioF() | reference point basis in ss3stk() (same logic, benchmark labels) |
# | SSplot_eqcurves() | ss3_eqcurve(), ss3_spm(), ss3plot_production() |
# | SSgetComps() | ss3_lencomp(), ss3_agecomp() |
# | (Pearson residuals in r4ss) | ss3_resid() (OSA), ss3plot_resid() |

# ---- About this guide ------------------------------------------------------------------------

# Part I (sections 3-9) goes through the diagnostics of an assessment, from the data to the
# prediction skill, with checks at each step (what to look at, what is a red flag, what to do
# about it), following the diagnostics cookbook. Part II (sections 10-14) covers stock status and
# its uncertainty, the timing of F and SSB in Kobe plots, the production function, season-, sex-
# and area-structured models (F by area, stock by area, partial impact of fleets) and the FLR
# outputs. Section 15 summarises the checks and the functions.

# The examples ship with the package:

# | Data | Model | Used for |
# | ss3simple | simulated two-sex model, one fishery and two surveys (r4ss "simple_small"; the example of the PIFSC ss3diags) | quick start, data, residuals |
# | ss3sma, retro.sma, aspm.sma | North Atlantic shortfin mako (cookbook) | runs tests, retrospective, status, Kobe |
# | ss3phk, retro.phk, aspm.phk | Pacific hake (cookbook) | Dirichlet-multinomial OSA residuals, hindcast |
# | pra.27.3a4a.rds | northern shrimp, four seasons, hermaphrodite, two areas, six area-specific fleets | structured models, F by area, impact |

####################################################################################################
# INSTALLATION
####################################################################################################

# ---- Required --------------------------------------------------------------------------------

# ss3diags2 needs FLCore, data.table and ggplot2. Install from r-universe (the same repositories
# provide the FLR packages):

if(FALSE){
  repos <- c(
    hw   = "https://henning-winker.r-universe.dev",
    flr  = "https://flr.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
  install.packages("ss3diags2", repos = repos)
}

# Or the development version from GitHub:

if(FALSE){
  remotes::install_github("Henning-Winker/ss3diags2")
}

# ---- Optional --------------------------------------------------------------------------------

# | Package | Needed for | Install |
# | r4ss | reading SS3 output (SS_output(), SSgetoutput(), SSsummarize()) | install.packages("r4ss") |
# | mvtnorm, MASS | MVLN draws (ss3mvln(), output = "iters") and Kobe contours | install.packages(c("mvtnorm", "MASS")) |
# | patchwork | the summary page ss3plot() and side-by-side plots | install.packages("patchwork") |
# | ss3diags | comparing with the original functions | remotes::install_github("jabbamodel/ss3diags") |
# | JABBA3, FLicc, FLRef | the models fed by the extracted inputs | see the packages |

# ---- Check the installation ------------------------------------------------------------------

library(ss3diags2)
library(ggplot2)
packageVersion("ss3diags2")
opt <- c("r4ss", "mvtnorm", "MASS", "patchwork", "ss3diags")
data.frame(package = opt, installed = vapply(opt, requireNamespace, logical(1), quietly = TRUE), row.names = NULL)

####################################################################################################
# QUICK START
####################################################################################################

# An SS3 run is read with r4ss::SS_output(); the report list is the input of ss3diags2. The
# simulated ss3simple model ships with the package (its SS3 files are in system.file("extdata",
# "simple_small", package = "r4ss")):

if(FALSE){
  ss3rep <- r4ss::SS_output("path/to/run", verbose = FALSE, printstats = FALSE)
}

data(ss3simple)
ss3_fleets(ss3simple)

# ---- Quick output: `ss3plot()` ---------------------------------------------------------------

# ss3plot() is the one-page summary of an assessment, as the JABBA3 summary plot:

# * top row: SSB (with the observed indices scaled to SSB as points), F and catch, with the
#   reference points (unfished in blue, targets in green, Bpa orange, Blim red);

# * middle row: SSB and F relative to their targets, and the Kobe plot with the uncertainty of
#   the terminal year;

# * bottom row: the stock-recruitment curve, the recruitment deviations with 95% intervals and
#   the production phase plot (catch against SSB on the Kobe colours, with the equilibrium yield
#   curve).

ss3plot(ss3simple, iters = 500)

# Uncertainty comes from MVLN draws of the SS3 Hessian. The targets are named after the benchmark
# SS3 used (here B40: 40% of unfished SSB, with Fmsy; section [Stock status]).

# ---- Stock status as FLR ---------------------------------------------------------------------

stk <- ss3stk(ss3simple)                     # FLStockR at the MLE, reference points in @refpts
stk@refpts
ss3plot_trj(stk)

####################################################################################################
# DATA: OBSERVATIONS AS FLR OBJECTS
####################################################################################################

# The extraction functions return the observations of an SS3 run, and with what = "exp" the model
# fit, as FLR objects. Each fleet is one element of an FLQuants; the SS3 structure goes into the
# FLQuant dimensions:

# | SS3 | FLQuant dimension | note |
# | length or age bin | len / age | data bins; empty bins of sampled compositions are filled with 0 |
# | year | year | |
# | sex (data Sexes 1-3) | unit | F, M; combined-sex data stay unique |
# | season | season | seasons stay seasons |
# | area | area | the area of the fleet |
# | partition 0/1/2 | fleet name | fleet.all, fleet.discard, fleet.retained |

ct  <- ss3_catch(ss3simple)                  # FLQuants: catch by fleet
idx <- ss3_indices(ss3simple)                # FLIndices of FLIndexBiomass: index, index.var (input SE)
lfd <- ss3_lencomp(ss3simple)                # FLQuants: length compositions by fleet, sex as unit
afd <- ss3_agecomp(ss3simple)

ss3plot_comps(lfd,type="annual",year=2015:2022)

# Every extractor also returns the long data.table (long = TRUE), with the FLQuant dimension
# columns:

head(ss3_lencomp(ss3simple, long = TRUE), 3)

# ---- Data overview ---------------------------------------------------------------------------

ss3plot_data(catch = ct, indices = idx, lencomp = lfd, agecomp = afd)

# ---- Compositions ----------------------------------------------------------------------------

# ss3plot_comps() shows compositions pooled over years, by year, as bubbles (only the years with
# data) or as mean length or age over time, split by unit, season or area (by = ), for a subset
# of fleets and units:

ss3plot_comps(afd, type = "bubble", by = "unit")
ss3plot_comps(lfd, type = "mean", by = "unit")

# ---- Inputs for JABBA3 and FLicc -------------------------------------------------------------

# ss2jabba3() returns catch in weight by fleet and annual indices with input SEs, the inputs of
# JABBA3; ss2flicc() returns length compositions by gear on common bins, with ss3_lhpar(), the
# inputs of FLicc.

j <- ss2jabba3(ss3simple)
names(j)

# Check 1: data. Look for contrast in the indices (a one-way trip gives little information on
# productivity), conflicting trends between indices, and compositions with few samples or gaps.
# ss3plot_nsamp() shows the input and effective sample sizes.

####################################################################################################
# LIFE HISTORY AND GROWTH
####################################################################################################

# ss3_lhpar() returns an FLPar of growth, maturity, length-weight and natural mortality, with a
# unit dimension when the model has two sexes; this is the lhpar of FLicc. Growth other than von
# Bertalanffy is approximated by a von Bertalanffy curve fitted to the SS3 growth at the end
# year.

lh <- ss3_lhpar(ss3simple)
lh
ss3plot_growth(ss3simple, lfd = lfd)

# Ridge plots show the length compositions by year, observed (filled) and expected (lines), with
# the growth curves by sex tracking cohorts through the bins:

ss3plot_ridges(lfd, fit = ss3_lencomp(ss3simple, what = "exp"), type = "both", fleets = "SURVEY1", lhpar = lh)

####################################################################################################
# CHECK 2: MODEL FITS
####################################################################################################

ss3plot_indices(idx, fit = ss3_indices(ss3simple, what = "exp"))

ss3plot_comps(lfd, fit = ss3_lencomp(ss3simple, what = "exp"), by = "unit")

# Check 2: fits. The fitted index should follow the trend of each index; systematic misfits at
# the start or end of a series, or an index fitted only where it agrees with another, point to
# conflicts between data sources. Pooled compositions should be fitted without systematic misfit
# of the small or large bins.

####################################################################################################
# CHECK 3: RESIDUAL DIAGNOSTICS
####################################################################################################

# ---- One-step-ahead residuals for compositions -----------------------------------------------

# Pearson residuals of multinomial compositions are correlated across bins (they must add up), so
# their patterns can mislead. One-step-ahead (OSA) residuals (Trijoulet et al. 2023) condition
# each bin on the bins before it; under a correct model they are independent N(0, 1). ss3_resid()
# computes them for the SS3 likelihood: multinomial with the adjusted sample size, or
# Dirichlet-multinomial with the precision implied by the SS3 effective sample size. Growth
# curves by sex (lhpar) track cohorts.

rl <- ss3_resid(ss3simple)                   # length; ss3_resid(ss3simple, "age") for ages
ss3plot_resid(rl, lhpar = lh)
ss3plot_resid(rl, type = "qq")

# Pacific hake uses the Dirichlet-multinomial for its age compositions:

data(ss3phk)
ss3plot_resid(ss3_resid(ss3phk, "age"), cohorts = FALSE)

# Reading OSA residuals. Bins are conditioned on the smaller bins: too few small fish gives
# negative residuals at small sizes and positive residuals at larger sizes. Diagonal patterns
# that follow the growth curves are year classes the model does not fit.

# ---- Runs tests and joint residuals ----------------------------------------------------------

# The runs test checks the log residuals of each index (or the mean length or age) for
# randomness; the 3-sigma limits come from the average moving range (Carvalho et al. 2021). North
# Atlantic shortfin mako:

data(ss3sma)
ss3_runs(ss3sma)
ss3plot_runs(ss3sma)

ss3plot_jres(ss3sma)

# Check 3: residuals. Red flags: a failed runs test (p < 0.05), residuals outside the 3-sigma
# limits, trends in the joint residuals (lowess line), RMSE above 30%. Actions: look for
# conflicts between indices, changes in catchability, selectivity or growth; never drop data only
# to pass a test.

####################################################################################################
# CHECK 4: RETROSPECTIVE ANALYSIS
####################################################################################################

# Retrospective runs refit the model with the last years of data removed (peels; r4ss::retro()),
# read with r4ss::SSgetoutput() (and optionally r4ss::SSsummarize()). ss3_retro() turns them into
# FLStocks (one FLStockR per peel); ss3_mohn() gives Mohn's rho and the forecast rho (the
# one-year-ahead projection of each peel against the reference run), as ss3diags::SShcbias().

data(retro.sma)
rt <- ss3_retro(retro.sma)
ss3_mohn(rt, c("SSB", "F"))

ss3plot_retro(rt, type = c("SSB", "F"))

# Check 4: retrospective. Mohn's rho outside -0.15 to 0.20 (long-lived species; -0.22 to 0.30 for
# short-lived species; Hurtado-Ferro et al. 2015) is a red flag, as are forecasts outside the
# confidence interval of the reference run. Actions: look at recent data, time-varying
# catchability or selectivity.

####################################################################################################
# CHECK 5: HINDCAST CROSS-VALIDATION
####################################################################################################

# The hindcast compares the observations of the peeled years with their one-year-ahead
# predictions from the retrospective runs (Kell et al. 2016). The mean absolute scaled error
# (MASE) scales the prediction error by that of a naive random-walk forecast: MASE < 1 means the
# model predicts better than the naive forecast. ss3_hindcast() does indices, mean length and
# mean age, as ss3diags::SSmase().

data(retro.phk)
hc <- ss3_hindcast(retro.phk)
hc$mase

ss3plot_hcxval(hc)

ss3plot_hcxval(ss3_hindcast(retro.sma))

# Check 5: prediction skill. MASE > 1 for an index means the model has no prediction skill for
# it. With indices that vary little the naive error is small and MAE.base.adj (default 0.1)
# limits the scaling.

####################################################################################################
# STOCK STATUS: FLSTOCKR AND MVLN
####################################################################################################

# ss3stk() builds a simplified FLStockR (as FLRef::ss2FLStockR()): SSB, F, recruitment and catch,
# with the reference points in @refpts. output = "mle" gives the maximum likelihood estimates
# with their standard errors, output = "iters" draws from a multivariate lognormal (MVLN)
# approximation (ss3mvln()):

# * method = "joint" (default) draws SSB, recruitment and F of all years and the reference points
#   from one MVLN with the full correlation matrix of the SS3 derived quantities (Kell et al. 2023;
#   the ss3om::ssmvln() approach). Each draw is a coherent trajectory and the reference points have
#   draws.

# * method = "annual" is the delta-MVLN of ss3diags::SSdeltaMVLN() (Walter and Winker 2020;
#   Winker et al. 2019), drawn year by year from the B ratio, F and the F reference.

stk <- ss3stk(ss3sma, output = "iters", iters = 1000, fwd = TRUE)
ss3plot_trj(stk)

# SSdeltaMVLN understates the B-F correlation. For the example models the year-by-year delta-MVLN
# gives a correlation between B/Bmsy and F/Fmsy in the terminal year of -0.24 for shortfin mako
# (SS3 covariance: -0.86) and -0.26 for Pacific hake (-0.76), so its Kobe clouds are rounder than
# the estimated uncertainty. The joint method keeps the estimated correlation; use it for advice.

# ---- Which targets: the reference point basis ------------------------------------------------

# SS3 can report stock status relative to different benchmarks, set in starter.ss (the B ratio
# denominator and the F report basis) and forecast.ss (the benchmarks; the MSY option can be set
# to the F at a biomass target or at an SPR target). The SS3 labels of the derived quantities
# also changed between versions. ss3stk() works out the basis with the logic of
# ss3diags::SSdeltaMVLN() and SSsettingsBratioF():

# 1. B ratio: the SS3 Bratio of the end year is compared with SSB relative to SSB virgin, SSB_MSY
# and SSB_Btgt (or SSB at F0.1); the closest defines the B target.

# 2. F basis: from F_std_basis (SS 3.30) or F_report_basis (3.24): absolute F, F/Fmsy, F/F at the
# B target, or F/F at the SPR target.

# 3. F reference: MSY, unless the MSY benchmark is the B target benchmark (then the F at the B
# target) or the SPR benchmark (then the SPR F); Fref = overrides it.

# The status values agree with SSdeltaMVLN() for all example models and for the SS3 benchmark
# options. The labels name the benchmark behind each target rather than the SS3 label, because
# "MSY" in SS3 output may be a proxy:

# | Basis | B target | F target | Catch reference |
# | MSY benchmark | Bmsy | Fmsy | MSY |
# | B ratio relative to B0, target fraction 40% | B40 | (from the F basis) | |
# | biomass target 40% (or MSY set to F(B40)) | B40 | F40 | Y40 (yield at B40) when MSY is the proxy |
# | SPR target 40% (or MSY set to F(SPR40)) | Bspr40 | Fspr40 | Yspr40 when MSY is the proxy |
# | F0.1 | BF0.1 | F0.1 | |

# The catch reference stays MSY when SS3 estimated MSY. The basis is stored in attr(stk,
# "basis"):

str(attr(ss3stk(ss3sma), "basis")[c("B", "F", "C", "bbasis", "fbasis")])
str(attr(ss3stk(ss3simple), "basis")[c("B", "F", "C", "bbasis", "fbasis")])

# ---- Key quantities as FLQuants --------------------------------------------------------------

q <- ss3_flqs(stk, c("SSB", "F", "B/Btgt", "F/Ftgt"))
names(q)

####################################################################################################
# KOBE PLOTS AND THE TIMING OF F AND SSB
####################################################################################################

# Critical for advice. SS3 reports SSB at spawning (typically at the start of the year), and F of
# year *y* acts during year *y*. The SSB that results from the fishing in year *y* is therefore
# SSB in year *y + 1*. The usual Kobe plot pairs F(*y*) with SSB(*y*): it shows the biomass the
# fishery in *y* started from, not the biomass it led to. When F is well above Fmsy, the stock
# status that follows from the current exploitation can be markedly worse than the plot of F(*y*)
# against SSB(*y*) suggests. Advice should state which convention is used and, when the two
# differ, show both.

# ss3plot_kobe(lag = 1) plots F(*y*) against B/Btgt(*y + 1*). The uncertainty cloud then pairs
# F(*Y*) with B(*Y + 1*) within each draw, which needs draws linked across years: use method =
# "joint". Build the stock with fwd = TRUE so that SSB at the start of the year after the last
# assessment year (an SS3 forecast quantity) is available; F is used up to the end year only.

k0 <- ss3plot_kobe(stk) + ggtitle("lag = 0: F(y) against B(y)")
k1 <- ss3plot_kobe(stk, lag = 1) + ggtitle("lag = 1: F(y) against B(y + 1)")
patchwork::wrap_plots(k0, k1)

# For North Atlantic shortfin mako the probability of being in the red quadrant (overfished and
# subject to overfishing) in 2015 is about 78% with F(2015) against SSB(2015) and about 93% with
# F(2015) against SSB(2016).

# Before using the lagged plot for advice. SSB(*y + 1*) for the last assessment year is a
# one-year forecast: it depends on the SS3 forecast settings (forecast.ss) and, where recruits
# mature within a year, on the forecast recruitment. Kell et al. (2023) pair F and SSB both ways
# in their figures; the convention is not standardised, which is why it must be stated with the
# advice.

####################################################################################################
# PRODUCTION FUNCTION AND EQUILIBRIUM
####################################################################################################

# ss3_eqcurve() extracts the SS3 equilibrium yield curve and ss3_spm() fits a Pella-Tomlinson
# production function to it, as a surplus production model equivalent of the age-structured model
# (the "production function" of JABBA3), and checks it against the SS3 benchmarks.

sp <- ss3_spm(ss3sma)
sp
ss3plot_production(ss3stk(ss3sma))
ss3plot_production(ss3stk(ss3sma), type = "dyn")

# Which unfished state? SS3 holds several: SSB_Virgin (start-year biology), SSB_unfished
# (benchmark years) and the dynamic B0 of each year. With time-varying weight-at-age or M
# (empirical wtatage, as for Pacific hake) these differ, and the equilibrium table can rest on
# yet another biology. ss3_spm() warns when its fit and the SS3 benchmarks differ by more than
# 5%.

####################################################################################################
# SEASON-, SEX- AND AREA-STRUCTURED MODELS
####################################################################################################

# The Stock Synthesis assessment of northern shrimp (*Pandalus borealis*) in Skagerrak and the
# Norwegian Deep (pra.27.3a4a, SS 3.30.22, 1970-2024) uses most of the structure SS3 offers: four
# seasons (spawning in season 1), two sexes with protandric hermaphroditism (born male, change to
# female), two areas (1: Skagerrak 3.a, 2: Norwegian Deep 4.a) with time-varying recruitment
# distribution, six area-specific fleets (Denmark, Norway, Sweden in each area) with retention
# and discards, and two surveys with sex-specific length compositions. F is reported as Fbar over
# ages 1-3 (Z - M) and the MSY benchmark is set to the F at B30%. The report and the SS3 input
# files ship with the package.

dir <- system.file("extdata", "pra.27.3a4a", package = "ss3diags2")
pra <- readRDS(file.path(dir, "pra.27.3a4a.rds"))
ss3_fleets(pra)

# ---- Observations by season, sex and area ----------------------------------------------------

# Because each fleet fishes in one area, the observations need no area split: the area is a
# property of the fleet.

ss3plot_catch(ss3_catch(pra), by = c("season", "area"))

lfp <- ss3_lencomp(pra)
ss3plot_comps(lfp, fit = ss3_lencomp(pra, what = "exp"), fleets = c("NORSURVEY3a", "NORSURVEY4a"), by = "unit")
ss3plot_resid(ss3_resid(pra, fleets = c("NORSURVEY3a", "NORSURVEY4a")), lhpar = ss3_lhpar(pra))

# ---- Selectivity, retention and the stock by area --------------------------------------------

ss3plot_sel(ss3_sel(pra), ret = ss3_sel(pra, "ret"))

ss3plot_byarea(ss3_byarea(pra))

# ---- Summary and status ----------------------------------------------------------------------

# The MSY benchmark is set to the F at B30%, so the targets are B30 and F30 and the catch
# reference is the yield at B30 (Y30):

ss3plot(pra, iters = 500, lag = 1)

# ---- F-at-age, partial F and F by area -------------------------------------------------------

# ss3_fa() returns F-at-age by fleet with the full structure (age, year, unit, season, area);
# F-at-age is on an annual scale within each season, and annual = TRUE sums the seasons weighted
# by their duration.

fa <- ss3_fa(pra)
ss3plot_fa(fa)

# Partial F splits the reported F (Fbar 1-3) among fleets, gears or countries by their catch in
# numbers over the F ages, so it adds up to F:

countries <- list(Denmark = c("DK3a", "DK4a"), Norway = c("NOR3a", "NOR4a"), Sweden = c("SWE3a", "SWE4a"))
ss3plot_partialF(ss3_partialF(pra, groups = countries))
ss3plot_partialF(ss3_partialF(pra), type = "line")

# Two notions of F by area. The partial F by area (groups = "area", method = "catch") is each
# area's share of the stock F and adds up to it. The F experienced by the fish in an area (method
# = "fbar", the Fbar of the fleets of that area) is the F of a stock by area; F's of different
# areas do not add up to the stock F when the stock is spread over the areas.

ss3plot_partialF(ss3_partialF(pra, groups = "area", method = "fbar"))

# ---- Stock by area ---------------------------------------------------------------------------

# ss3stk(area = TRUE) returns FLStocks with the stock and one FLStockR per area: SSB, recruitment
# and catch by area, and the F experienced in each area. Area reference points are
# approximations: B0 and R0 are the unfished values of the area, the B target is scaled by the
# area's share of B0 and the F target is that of the stock.

sa <- ss3stk(pra, area = TRUE)
ss3plot_trj(sa)

# ---- Partial impact of fleets ----------------------------------------------------------------

# ss3_impact() follows Wang et al. (2009): SSB is projected without fishing and with the fleet
# groups added one by one, so the impact of a group is the SSB it removes, and all impacts add up
# to the dynamic B0 minus SSB. It runs on the SS3 report without re-running SS3: annual total
# mortality from the SS3 numbers at age (seasons, areas and sexes pooled), fishing mortality $Z -
# M$ split among fleets by their catch at age, the SS3 spawning output at age and recruitment
# deviations, and each scenario starting from the initial equilibrium under the initial F of its
# fleets. The run with all fleets reproduces the SS3 SSB and the run without fishing the SS3
# dynamic B0:

imp <- ss3_impact(pra, groups = countries)
imp
ss3plot_impact(imp)
ss3plot_impact(imp, type = "prop")

# Approximations. When fleets are removed, movement, sex change and the timing within the year
# stay at their fitted values. The order in which the groups are added matters for the share of
# each group when the stock-recruit relationship is not flat; type = "remove" gives the loss of
# SSB caused by each group with all others fishing.

####################################################################################################
# FLR OUTPUTS
####################################################################################################

# | Output | Function | Class |
# | catch, discards by fleet | ss3_catch(), ss3_discards() | FLQuants |
# | indices | ss3_indices() | FLIndices (FLIndexBiomass) |
# | length, age compositions | ss3_lencomp(), ss3_agecomp() | FLQuants |
# | OSA / Pearson residuals | ss3_resid() | FLQuants |
# | life history | ss3_lhpar() | FLPar (unit by sex) |
# | selectivity, retention | ss3_sel() | FLQuants |
# | F-at-age by fleet | ss3_fa() | FLQuants (age, year, unit, season, area) |
# | partial F | ss3_partialF() | FLQuants (total F as attribute) |
# | stock by area | ss3_byarea(); ss3stk(area = TRUE) | FLQuants; FLStocks |
# | stock status | ss3stk() | FLStockR (@refpts) |
# | retrospective peels | ss3_retro() | FLStocks |
# | key quantities | ss3_flqs() | FLQuants |
# | MVLN draws | ss3mvln()$flqs, ss3_flqs(mvln) | FLQuants with iterations |

# long2flq(), long2flqs() and flqs2long() convert between the long tables and FLR objects. The
# stock objects work with FLRef and ggplotFL:

library(ggplotFL)
plot(ss3_flqs(ss3stk(ss3sma), c("SSB", "F", "Catch")))

# ---- MVLN draws as FLQuants with iterations --------------------------------------------------

# ss3mvln() returns the draws as FLQuants with one iteration per draw ($flqs): SSB, F,
# recruitment, catch, B/Btgt and F/Ftgt, and B/B0 and Catch/MSY where the reference points are
# reported. ss3_flqs() selects quantities from it, and ss3plot_flqs() plots the median with 50%
# and 90% bands, with the ratio panels named after the benchmark and a few single draws ("worms")
# if wanted:

mv <- ss3mvln(ss3sma, mc = 1000)              # joint MVLN draws
mv$flqs
dim(mv$flqs$SSB)                              # 1000 iterations
ss3plot_flqs(ss3_flqs(mv, c("SSB", "F", "B/Btgt", "F/Ftgt")), worms = 5)

# The year-by-year delta-MVLN draws (method = "annual") as FLQuants, for Pandalus with the B30
# and F30 targets:

mp <- ss3mvln(pra, mc = 1000, method = "annual")
ss3plot_flqs(ss3_flqs(mp, c("SSB", "Rec", "B/Btgt", "F/Ftgt")))

# The same FLQuants come from the stock with iterations, ss3_flqs(ss3stk(ss3sma, output =
# "iters")), and work with ggplotFL (plot(mv$flqs)).

####################################################################################################
# SUMMARY OF CHECKS AND FUNCTIONS
####################################################################################################

# | # | Check | Red flag | Action |
# | 1 | data (ss3plot_data, ss3plot_nsamp) | no contrast; conflicting indices; few samples | sensitivity runs; data weighting |
# | 2 | fits (ss3plot_indices, ss3plot_comps) | systematic misfit | conflicts between data; selectivity, catchability |
# | 3 | residuals (ss3plot_resid, ss3plot_runs, ss3plot_jres) | failed runs test; residuals outside 3 sigma; RMSE > 30% | time-varying q, selectivity or growth |
# | 4 | retrospective (ss3_mohn, ss3plot_retro) | Mohn's rho outside -0.15 to 0.20 (short-lived: -0.22 to 0.30) | recent data; time-varying processes |
# | 5 | hindcast (ss3_hindcast, ss3plot_hcxval) | MASE > 1 | index without prediction skill |
# | 6 | status (ss3stk, ss3plot_kobe) | Kobe lag changes the quadrant | state the convention; show both |

# ---- Function reference ----------------------------------------------------------------------

# | Area | Functions |
# | Example data | ss3simple, ss3sma, retro.sma, aspm.sma, ss3phk, retro.phk, aspm.phk, extdata/pra.27.3a4a |
# | Observations | ss3_fleets(), ss3_catch(), ss3_discards(), ss3_indices(), ss3_lencomp(), ss3_agecomp() |
# | Model inputs | ss2jabba3(), ss2flicc(), ss3_lhpar() |
# | Observation plots | ss3plot_data(), ss3plot_catch(), ss3plot_indices(), ss3plot_comps(), ss3plot_nsamp(), ss3plot_growth(), ss3plot_ridges() |
# | Residuals | ss3_resid(), ss3plot_resid(), ss3_runs(), ss3plot_runs(), ss3plot_jres(), ss3_iresid() |
# | Retro and hindcast | ss3_retro(), ss3_mohn(), ss3plot_retro(), ss3_hindcast(), ss3plot_hcxval() |
# | Status | ss3stk(), ss3mvln(), ss3_flqs(), ss3plot_flqs(), ss3plot(), ss3plot_trj(), ss3plot_kobe() |
# | Production | ss3_eqcurve(), ss3_spm(), ss3plot_production() |
# | Structure | ss3_sel(), ss3plot_sel(), ss3_byarea(), ss3plot_byarea(), ss3plot_sr() |
# | F and impact | ss3_fa(), ss3plot_fa(), ss3_partialF(), ss3plot_partialF(), ss3_impact(), ss3plot_impact() |
# | Utilities | long2flq(), long2flqs(), flqs2long(), theme_ss3(), ss3cols() |

