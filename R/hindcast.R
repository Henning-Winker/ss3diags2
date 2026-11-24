# hindcast.R - hindcast cross-validation (one-year-ahead predictions of indices, mean length / age) and MASE

#' Observations and fits of one run for the hindcast
#'
#' Observations and fits of one run: index (cpue) or mean length / age by fleet and year
#'
#' @param ss3rep SS3 report (\code{r4ss::SS_output()})
#' @param quant data type: \code{"cpue"}, \code{"len"}, \code{"age"}, \code{"con"} (conditional age-at-length) or \code{"size"}
#' @return data.table of index (cpue) or mean length / age by fleet and year
#' @noRd
.hc_obs <- function(ss3rep, quant){
  if(quant == "cpue"){
    cp <- as.data.table(ss3rep$cpue)
    d <- data.table(fleet = cp$Fleet, name = .fleetnames(ss3rep)[cp$Fleet], year = cp$Yr, season = cp$Seas,
                    obs = cp$Obs, exp = cp$Exp, se = cp$SE, use = .col(cp, "Use", 1), like = cp$Like)
    return(d)
  }
  db <- .comp_db(ss3rep, quant, ghost = TRUE)
  if(quant == "con"){
    # conditional age-at-length: mean age by length bin, combined over length bins weighted by N (as SScompsTA1.8)
    lb <- db[, .(o = sum(Bin*Obs)/sum(Obs), e = sum(Bin*Exp)/sum(Exp), n = mean(Nsamp_adj),
                 like = if(all(is.na(Like))) NA_real_ else 1), by = .(Fleet, Yr, Seas, Lbin_lo)]
    d <- lb[is.finite(o) & is.finite(e), .(obs = sum(o*n)/sum(n), exp = sum(e*n)/sum(n),
                                          like = if(all(is.na(like))) NA_real_ else 1), by = .(Fleet, Yr, Seas)]
  } else {
    # mean length / age / size of each fleet x year x season (sexes and partitions pooled), as SScompsTA1.8, with
    # the SE of the observed mean from the expected variance and the (adjusted or Dirichlet-multinomial) effective
    # sample size (Francis 2011, method TA1.8)
    if(!"DM_effN" %in% names(db)) db[, DM_effN := NA_real_]
    db[, neff := ifelse(is.finite(DM_effN), DM_effN, Nsamp_adj)]
    d <- db[, .(obs = sum(Bin*Obs)/sum(Obs), exp = sum(Bin*Exp)/sum(Exp),
                v = sum(Bin^2*Exp)/sum(Exp) - (sum(Bin*Exp)/sum(Exp))^2,
                n = sum(unique(data.table(Sex = .col(.SD, "Sex", 1), Part, neff))$neff),
                like = if(all(is.na(Like))) NA_real_ else 1), by = .(Fleet, Yr, Seas)]
    d[, sd := sqrt(pmax(v, 0)/n)]
  }
  data.table(fleet = d$Fleet, name = .fleetnames(ss3rep)[d$Fleet], year = d$Yr, season = d$Seas, obs = d$obs,
             exp = d$exp, se = if(quant == "con") NA_real_ else d$sd/d$obs, use = 1, like = d$like)
}

# composition data base: len, age, con (conditional age-at-length) or size; ghost = TRUE adds the excluded
#' Composition database of a data type, optionally with ghost observations
#'
#' (ghost) observations, e.g. the peeled years of retrospective runs
#'
#' @param ss3rep SS3 report (\code{r4ss::SS_output()})
#' @param quant data type: \code{"cpue"}, \code{"len"}, \code{"age"}, \code{"con"} (conditional age-at-length) or \code{"size"}
#' @param ghost include ghost observations (e.g. the peeled years of retrospective runs)
#' @return data.table
#' @noRd
.comp_db <- function(ss3rep, quant, ghost = FALSE){
  nm <- c(len = "lendbase", age = "agedbase", con = "condbase", size = "sizedbase")[[quant]]
  db <- as.data.table(ss3rep[[nm]])
  if(ghost){
    g <- as.data.table(ss3rep[[paste0("ghost", nm)]])
    if(nrow(g)){ g[, Fleet := abs(Fleet)]; db <- rbind(db, g, fill = TRUE) }
  }
  if(!nrow(db)) stop("no ", c(len = "length", age = "age", con = "conditional age-at-length", size = "generalised size")[[quant]],
                     " compositions in this SS3 report", call. = FALSE)
  db
}

#' Hindcast cross-validation
#'
#' One-year-ahead hindcasts of the abundance indices (or mean length / age) from retrospective runs: each peel with end
#' year \eqn{T} predicts the observation in \eqn{T + 1}; the prediction skill is the mean absolute scaled error, MASE
#' = MAE of the log prediction residuals / MAE of a naive (random walk) prediction from the previous evaluation year,
#' as \code{ss3diags::SSmase()} (Kell et al. 2021, Front. Mar. Sci. 8: 644563). MASE < 1: the model predicts better
#' than the naive forecast. \code{MASE.adj} uses \code{max(MAE.naive, MAE.base.adj)} as the scaler, to avoid inflated
#' MASE for nearly flat series.
#'
#' @param retro list of SS3 reports (reference run first, e.g. \code{r4ss::SSgetoutput()}), or their
#'   \code{r4ss::SSsummarize()} summary (indices only for the hindcast)
#' @param peels years peeled from each run (default \code{0:(n - 1)}); needed only when the reports do not give the
#'   peel end year (\code{Retro_year}), e.g. for a summary
#' @param quant \code{"cpue"} (default), \code{"len"}, \code{"age"}, \code{"con"} (conditional age-at-length:
#'   mean age by length bin combined over bins) or \code{"size"} (generalised size compositions); for compositions
#'   the observations excluded from the peels (SS3 ghost data) are included
#' @param MAE.base.adj minimum naive MAE for \code{MASE.adj}
#' @param fleets optional fleets (numbers or names)
#' @return list of class \code{ss3hc}: \code{data} (data.table: name, year, obs, se, ref (reference fit), peel,
#'   exp (peel fit), endyr, pred (TRUE for the one-year-ahead prediction)), \code{mase} (data.table: index, MASE,
#'   MAE.PR, MAE.base, MASE.adj, n.eval, plus "joint")
#' @examples
#' data(retro.sma)
#' hc <- ss3_hindcast(retro.sma)
#' hc$mase
#' ss3plot_hcxval(hc)
#' @export
ss3_hindcast <- function(retro, quant = c("cpue", "len", "age", "con", "size"), MAE.base.adj = 0.1, fleets = NULL,
                         peels = NULL){
  quant <- match.arg(quant)
  if(.is_summary(retro) && quant != "cpue")
    stop("SSsummarize() output holds the indices only: quant = \"cpue\" (full reports for compositions)",
         call. = FALSE)
  retro <- .retro_reps(retro, peels)
  ey <- vapply(retro, .endyr, 0)
  runs <- rbindlist(lapply(seq_along(retro), function(i) .hc_obs(retro[[i]], quant)[, `:=`(peel = as.character(ey[i]),
                                                                                         endyr = ey[i])]))
  ref <- runs[peel == as.character(ey[1])]
  ref <- ref[use > 0 & !is.na(like)]
  if(!is.null(fleets)){ fn <- if(is.numeric(fleets)) .fleetnames(retro[[1]])[fleets] else fleets; ref <- ref[name %in% fn] }
  # one observation per year: all of them, whatever the season (e.g. a survey whose month varies between years);
  # several per year (seasonal index): the first season only, as ss3diags::SSmase(Season = "default")
  ref <- ref[, if(anyDuplicated(year)) .SD[season == season[1]] else .SD, by = name]
  eyrs <- sort(ey)
  dat <- rbindlist(lapply(unique(ref$name), function(f){
    r <- ref[name == f]
    p <- runs[name == f & paste(year, season) %in% paste(r$year, r$season) & peel != as.character(ey[1]) &
                year <= endyr + 1]
    p[, pred := year == endyr + 1]
    rbind(r[, .(name, year, obs, se, ref = exp, peel = "Ref", exp = exp, endyr = ey[1], pred = FALSE)],
          merge(p[, .(name, year, peel, exp, endyr, pred)], r[, .(year, obs, se, ref = exp)], by = "year"),
          use.names = TRUE)
  }))
  mase <- rbindlist(lapply(unique(ref$name), function(f){
    r <- ref[name == f]
    have <- eyrs[eyrs %in% r$year]
    pe <- have[-1]
    pr <- dat[name == f & pred & year %in% pe]
    if(length(pe) < 1 || !nrow(pr))
      return(data.table(index = f, MASE = NA_real_, MAE.PR = NA_real_, MAE.base = NA_real_, MASE.adj = NA_real_, n.eval = 0L))
    o <- r$obs[match(have, r$year)]
    naive <- log(o[-1]) - log(o[-length(o)])
    maepr <- mean(abs(log(pr$exp) - log(pr$obs))); scaler <- mean(abs(naive))
    data.table(index = f, MASE = maepr/scaler, MAE.PR = maepr, MAE.base = scaler,
               MASE.adj = maepr/max(scaler, MAE.base.adj), n.eval = nrow(pr))
  }))
  pr <- dat[pred == TRUE]
  ok <- mase[is.finite(MASE)]
  if(nrow(ok)){
    # joint: all prediction and naive residuals pooled
    j <- rbindlist(lapply(ok$index, function(f){
      r <- ref[name == f]; have <- eyrs[eyrs %in% r$year]; o <- r$obs[match(have, r$year)]
      list(pr = abs(log(pr[name == f & year %in% have[-1]]$exp) - log(pr[name == f & year %in% have[-1]]$obs)),
           nv = abs(log(o[-1]) - log(o[-length(o)])))
    }))
    jp <- mean(unlist(j$pr)); jn <- mean(unlist(j$nv))
    mase <- rbind(mase, data.table(index = "joint", MASE = jp/jn, MAE.PR = jp, MAE.base = jn,
                                   MASE.adj = jp/max(jn, MAE.base.adj), n.eval = sum(ok$n.eval)))
  }
  structure(list(data = dat[], mase = mase[], quant = quant, endyrs = ey), class = "ss3hc")
}

#' @export
print.ss3hc <- function(x, ...){
  cat("Hindcast cross-validation (", x$quant, "), ", length(x$endyrs) - 1, " peels. MASE:\n", sep = "")
  print(x$mase, digits = 3)
  invisible(x)
}

#' Plot hindcast cross-validation
#'
#' Observations (with lognormal intervals for indices), the reference fit, the peel fits to their end year with the
#' one-year-ahead prediction as a dashed segment ending in a point, and the predicted observations coloured by peel;
#' MASE in the panel titles (MASE.adj in brackets when smaller), as \code{ss3diags::SSplotHCxval()} and
#' \code{JABBA3::jbplot_hcxval()}.
#'
#' @param x output of \code{\link{ss3_hindcast}()}, or a list of SS3 reports
#' @param ci confidence level of the observation intervals: indices from their log-scale SE; mean length / age
#'   from the SE of the mean (expected variance / effective sample size, Francis 2011 method TA1.8)
#' @param xmin first year shown (default 5 years before the first peel end year)
#' @param ncol facet columns
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' data(retro.sma)
#' ss3plot_hcxval(ss3_hindcast(retro.sma), xmin = 2000)
#' @export
ss3plot_hcxval <- function(x, ci = 0.95, xmin = NULL, ncol = NULL, text = TRUE){
  .drop_text((function(){
  if(!inherits(x, "ss3hc")) x <- ss3_hindcast(x)
  d <- copy(x$data); m <- x$mase[index != "joint"]
  lab <- m[, .(index, lab = ifelse(is.finite(MASE), paste0(index, ": MASE = ", sprintf("%.2f", MASE),
                                   ifelse(MASE.adj < MASE, sprintf(" (%.2f)", MASE.adj), "")),
                                   paste0(index, ": no prediction")))]
  d[, panel := factor(lab$lab[match(name, lab$index)], levels = lab$lab)]
  xmin <- xmin %||% (min(x$endyrs) - 5)
  d <- d[year >= xmin]
  pl <- c("Ref", as.character(x$endyrs[-1]))
  cols <- stats::setNames(c("black", ss3cols(length(pl) - 1)), pl)
  o <- unique(d[peel == "Ref", .(panel, year, obs, se)])   # intervals of the reference run (all observations)
  z <- stats::qnorm(1 - (1 - ci)/2)
  if(x$quant == "cpue") o[, `:=`(lo = obs*exp(-z*se), hi = obs*exp(z*se))] else
    o[, `:=`(lo = obs*(1 - z*se), hi = obs*(1 + z*se))]   # normal: se is the CV of the observed mean
  pe <- d[peel != "Ref"][order(panel, peel, year)]
  fit <- pe[year <= endyr]; fc <- pe[year == endyr | year == endyr + 1][, if(.N == 2) .SD, by = .(panel, peel)]
  pp <- d[pred == TRUE]
  p <- ggplot2::ggplot(mapping = ggplot2::aes(year))
  if(any(is.finite(o$lo))) p <- p + ggplot2::geom_ribbon(data = o[is.finite(lo)], ggplot2::aes(ymin = lo, ymax = hi),
                                                         fill = "grey80", alpha = 0.6)
  p + ggplot2::geom_point(data = o, ggplot2::aes(y = obs), shape = 21, fill = "white", size = 2) +
    ggplot2::geom_line(data = d[peel == "Ref"], ggplot2::aes(y = exp, colour = peel), linewidth = 0.8) +
    ggplot2::geom_line(data = fit, ggplot2::aes(y = exp, colour = peel), linewidth = 0.5) +
    ggplot2::geom_line(data = fc, ggplot2::aes(y = exp, colour = peel), linetype = 2, linewidth = 0.5) +
    ggplot2::geom_point(data = pp, ggplot2::aes(y = exp, fill = peel), shape = 21, size = 1.8, colour = "black") +
    ggplot2::geom_point(data = pp, ggplot2::aes(y = obs, fill = peel), shape = 21, size = 2.6, colour = "black") +
    ggplot2::scale_colour_manual(values = cols, breaks = pl) + ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = ncol) + ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "Year", y = switch(x$quant, cpue = "Index", len = "Mean length", age = "Mean age",
                                                 con = "Mean age (conditional on length)", size = "Mean size")) +
    scale_x_year() + theme_ss3()
  })(), text)
}

utils::globalVariables(c("Bin", "Obs", "Exp", "Like", "Nsamp_adj", "like", "pred", "index", "MASE", "MASE.adj",
                         "season", "lo", "hi", "DM_effN", "neff", "Part", "v", "sd"))
