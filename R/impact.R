# impact.R - partial impact of fleets on the spawning stock (Wang et al. 2009)

# numbers at age at the beginning of season s, summed over areas, sexes and morphs: matrix ages x years
.natage_mat <- function(ss3rep, era = c("TIME", "FORE"), seas = 1){
  na <- as.data.table(ss3rep$natage)
  bm <- .col(na, c("Beg/Mid", "BegMid"))
  na <- na[bm == "B" & Era %in% era & Seas == seas]
  ages <- names(na)[grepl("^[0-9]+$", names(na))]
  d <- na[, lapply(.SD, sum), by = Yr, .SDcols = ages]
  m <- t(as.matrix(d[, ages, with = FALSE])); colnames(m) <- d$Yr
  m[, order(d$Yr), drop = FALSE]
}

# spawning output at age (females at spawning time x Mat*Fecund), summed over areas: matrix ages x years
.spawn_contrib <- function(ss3rep, years){
  na <- as.data.table(ss3rep$natage); bm <- .col(na, c("Beg/Mid", "BegMid"))
  sp <- ss3rep$spawnseas %||% 1
  ages <- names(na)[grepl("^[0-9]+$", names(na))]
  na <- na[bm == "B" & Era == "TIME" & Seas == sp & Yr %in% years & (Sex == 1 | (ss3rep$nsexes %||% 1) == 1)]
  eg <- as.data.table(ss3rep$endgrowth)[Seas == sp]
  mf <- .col(eg, c("Mat*Fecund", "Mat_F_wtatage"))
  fec <- data.table(fec = mf, Morph = eg$Morph, age = as.character(.col(eg, c("int_Age", "Age_Beg"))))
  d <- melt(na, id.vars = c("Yr", "Morph"), measure.vars = ages, variable.name = "age", variable.factor = FALSE)
  d <- merge(d, unique(fec[, .(Morph, age, fec)], by = c("Morph", "age")), by = c("Morph", "age"), all.x = TRUE)
  d <- d[, .(s = sum(value*fec, na.rm = TRUE)), by = .(Yr, age)]
  m <- as.matrix(dcast(d, age ~ Yr, value.var = "s")[order(as.numeric(age))][, -1])
  m[, as.character(years), drop = FALSE]
}

# SS3 annual natural mortality at age (seasons weighted by duration, sexes averaged); NULL if not reported
.ss3_M <- function(ss3rep, ages){
  eg <- as.data.table(ss3rep$endgrowth)
  if(!nrow(eg) || !"M" %in% names(eg)) return(NULL)
  dur <- ss3rep$seasdurations %||% 1
  eg[, age := .col(eg, c("int_Age", "Age_Beg"))]
  m <- eg[, .(M = sum(M*dur[Seas])), by = .(Morph, age)][, .(M = mean(M)), by = age]
  m$M[match(as.numeric(ages), m$age)]
}

# stock-recruit shape f(S) (R0 cancels in the ratio f(S')/f(S))
.sr_shape <- function(ss3rep, sr){
  if(sr == "constant") return(function(S) rep(1, length(S)))
  p <- as.data.table(ss3rep$parameters)
  h <- p[Label %in% c("SR_BH_steep", "SR_Ricker_beta", "SR_steep", "SR_BH_flat_steep")]$Value[1]
  S0 <- ss3rep$SBzero %||% NA
  type <- ss3rep$SRRtype %||% 3
  if(is.na(h) || is.na(S0)) { warning("SR parameters not found: constant recruitment assumed", call. = FALSE)
    return(function(S) rep(1, length(S))) }
  switch(as.character(type),
         "3" = , "6" = function(S) S/(S0*(1 - h) + S*(5*h - 1)),
         "2" = function(S) S*exp(-h*S/S0),
         { warning("SR type ", type, " not supported: constant recruitment assumed", call. = FALSE)
           function(S) rep(1, length(S)) })
}

#' Partial impact of fleets on the spawning stock
#'
#' Computes the impact of fleets or fleet groups on the spawning stock after Wang et al. (2009): the stock is
#' projected through the assessment years without fishing and with the groups added one at a time
#' (\code{type = "cumulative"}), and the impact of a group is the loss of SSB when it is added. The SSB with no
#' fishing is the dynamic unfished SSB (dynamic B0) and the SSB with all fleets is the SS3 estimate, so the impacts
#' of all groups add up to dynamic B0 minus SSB.
#'
#' The projection works on the SS3 report and needs no re-run of SS3. Total mortality at age per year is taken from
#' the SS3 numbers at age (summed over seasons, areas and sexes: \eqn{Z = \log(N_{a,y}/N_{a+1,y+1})}) and split into
#' fishing mortality \eqn{Z - M} (with the SS3 natural mortality M) is split among fleets in proportion to their catch
#' at age in numbers. This keeps the seasonal, spatial and sex structure of SS3 in the realised annual rates, so the
#' projection with all fleets reproduces the SS3 numbers at age and SSB. Recruitment follows the SS3 stock-recruit
#' curve with the estimated deviations kept (\eqn{R' = R\,f(S')/f(S)}), and spawning output at age is that of SS3
#' (females, maturity x fecundity, at spawning time). Removing fleets changes F, so movement, sex change and
#' within-year timing are held at their fitted values: an approximation that is close for moderate changes in F.
#' With \code{init = "virgin"}, the no-fishing projection reproduces the SS3 dynamic B0 (\code{ss3rep$Dynamic_Bzero},
#' reported in \code{$check}).
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param groups \code{NULL} (each catch fleet), \code{"area"} or a named list of fleet names or numbers (see
#'   \code{\link{ss3_partialF}()}); the order of the list is the order in which groups are added
#' @param type \code{"cumulative"} (default, Wang et al. 2009: groups added one at a time, impacts add up) or
#'   \code{"remove"} (each group removed from the fishery in turn; impacts do not add up)
#' @param init initial numbers: \code{"virgin"} (default; unfished numbers at age in the first year, as the SS3
#'   dynamic B0, which the no-fishing scenario then reproduces) or \code{"model"} (the SS3 numbers in the first year,
#'   which include the initial depletion). With \code{"virgin"} the scenario with all fleets starts unfished and joins
#'   the SS3 SSB after about a generation; with \code{"model"} it equals the SS3 SSB from the first year but the
#'   impact of fishing before the first year is not counted.
#' @param sr recruitment: \code{"ss3"} (SS3 stock-recruit curve with deviations, default) or \code{"constant"} (the
#'   SS3 recruitments)
#' @param M natural mortality: \code{"ss3"} (default; SS3 M at age, seasons weighted by duration, with the fishing
#'   mortality \eqn{Z - M} split among fleets by their catch at age in numbers) or \code{"baranov"} (M is the
#'   remainder of \eqn{Z} after \eqn{F_f = Z\,C_f/D}; needs no biology but is biased when F is seasonal)
#' @param years years to project (default: start to end year)
#' @return object of class \code{ss3impact}: \code{ssb} (data.table year, scenario, ssb, rec), \code{impact}
#'   (year, group, impact, prop), \code{flqs} (FLQuants of SSB by scenario), \code{check} (SS3 versus projected SSB
#'   with all fleets, and the SS3 dynamic B0 when reported) and \code{groups}
#' @references Wang, S.-P., Maunder, M.N., Aires-da-Silva, A., Bayliff, W.H. (2009). Evaluating fishery impacts:
#'   application to bigeye tuna (Thunnus obesus) in the eastern Pacific Ocean. Fisheries Research 99: 106-111.
#' @seealso \code{\link{ss3plot_impact}()}, \code{\link{ss3_partialF}()}
#' @examples
#' \dontrun{
#' imp <- ss3_impact(ss3rep, groups = list(Denmark = 1:2, Norway = 3:4, Sweden = 5:6))
#' ss3plot_impact(imp)
#' }
#' @export
ss3_impact <- function(ss3rep, groups = NULL, type = c("cumulative", "remove"), init = c("virgin", "model"),
                       sr = c("ss3", "constant"), M = c("ss3", "baranov"), years = NULL){
  type <- match.arg(type); init <- match.arg(init); sr <- match.arg(sr); M <- match.arg(M)
  fl <- ss3_fleets(ss3rep)
  fleets <- fl[type %in% c("catch", "bycatch")]$name
  grp <- .fleet_groups(ss3rep, groups, fleets)
  yrs <- years %||% seq(ss3rep$startyr, ss3rep$endyr)
  N <- .natage_mat(ss3rep)
  if(!all(as.character(c(yrs, max(yrs) + 1)) %in% colnames(N))){
    yrs <- yrs[as.character(yrs + 1) %in% colnames(N)]
    if(!length(yrs)) stop("numbers at age for the year after the end year not found", call. = FALSE)
  }
  A <- nrow(N); ny <- length(yrs); ya <- as.character(yrs)
  N1 <- N[, as.character(yrs + 1), drop = FALSE]; N <- N[, ya, drop = FALSE]
  # total mortality: Z and deaths by age, plus group pooled with the age before
  nxt <- rbind(N1[-1, , drop = FALSE], 0)
  Z <- log(N/nxt); D <- N - nxt
  pg <- c(A - 1, A)
  Z[pg, ] <- rep(log(colSums(N[pg, , drop = FALSE])/N1[A, ]), each = 2)
  D[pg, ] <- rep(colSums(N[pg, , drop = FALSE]) - N1[A, ], each = 2)
  # fleet F at age from the catch-at-age shares of the deaths
  ca <- .catage_annual(ss3rep)[year %in% yrs]
  nm <- .fleetnames(ss3rep)
  Ff <- lapply(stats::setNames(fleets, fleets), function(f){
    m <- matrix(0, A, ny, dimnames = dimnames(N))
    x <- ca[fleet == match(f, nm)]
    m[cbind(match(x$age, as.numeric(rownames(N))), match(x$year, yrs))] <- x$data
    m[pg, ] <- rep(colSums(m[pg, , drop = FALSE]), each = 2)
    m
  })
  Ctot <- Reduce(`+`, Ff)
  bad <- !is.finite(Z) | Z <= 0 | D <= 0
  Mb <- if(M == "ss3") .ss3_M(ss3rep, rownames(N)) else NULL
  if(is.null(Mb)){
    # Baranov: F_f = Z C_f / D, M the remainder
    Ff <- lapply(Ff, function(m){ r <- m/D*Z; r[bad | !is.finite(r)] <- 0; r })
    Mm <- Z - Reduce(`+`, Ff)
    if(any(Mm[!bad] < -1e-6)) warning("catch exceeds deaths for some ages/years: M set to >= 0", call. = FALSE)
    Mm[bad] <- 0; Mm <- pmax(Mm, 0)
  } else {
    # SS3 natural mortality; fishing mortality Z - M split among fleets by catch at age in numbers
    Mm <- matrix(Mb, A, ny)
    Ft <- pmax(Z - Mm, 0); Ft[bad] <- 0
    Ff <- lapply(Ff, function(m){ r <- Ft*m/Ctot; r[!is.finite(r)] <- 0; r })
    Mm <- Z - Reduce(`+`, Ff); Mm[bad] <- matrix(Mb, A, ny)[bad]
  }
  Fall <- Reduce(`+`, Ff)
  # spawning output per fish at the beginning of the year, scaled to the SS3 SSB
  fec <- .spawn_contrib(ss3rep, yrs)/N
  fec[!is.finite(fec)] <- 0
  ts <- as.data.table(ss3rep$timeseries)
  ssb3 <- ts[Yr %in% yrs & !is.na(SpawnBio), .(s = sum(SpawnBio)), by = Yr][match(yrs, Yr)]$s
  fec <- sweep(fec, 2, ssb3/colSums(fec*N), "*")
  dur <- ss3rep$seasdurations %||% 1
  tsp <- sum(dur[seq_len((ss3rep$spawnseas %||% 1) - 1)])
  fS <- .sr_shape(ss3rep, sr)
  Nv <- if(init == "virgin") { v <- .natage_mat(ss3rep, era = "VIRG"); v[, ncol(v)] } else N[, 1]
  proj <- function(Fsum){
    Np <- matrix(0, A, ny); Np[, 1] <- Nv; S <- R <- numeric(ny)
    for(i in seq_len(ny)){
      Zp <- Mm[, i] + Fsum[, i]
      S[i] <- sum(Np[, i]*fec[, i]*exp(-(Zp - Z[, i])*tsp))
      R[i] <- N[1, i]*fS(S[i])/fS(ssb3[i])
      Np[1, i] <- R[i]
      if(i < ny){
        s <- Np[, i]*exp(-Zp)
        Np[-1, i + 1] <- s[-A]; Np[A, i + 1] <- Np[A, i + 1] + s[A]
      }
    }
    list(ssb = S, rec = R)
  }
  zero <- 0*Fall
  gsum <- lapply(grp, function(g) Reduce(`+`, Ff[g]))
  sc <- if(type == "cumulative"){
    lv <- c("no fishing", vapply(seq_along(grp), function(k) paste0("+", names(grp)[k]), ""))
    stats::setNames(lapply(0:length(grp), function(k) if(k == 0) zero else Reduce(`+`, gsum[seq_len(k)])), lv)
  } else {
    lv <- c("no fishing", paste0("-", names(grp)), "all fleets")
    stats::setNames(c(list(zero), lapply(names(grp), function(g) Fall - gsum[[g]]), list(Fall)), lv)
  }
  res <- rbindlist(lapply(names(sc), function(s){ p <- proj(sc[[s]]); data.table(year = yrs, scenario = s,
                                                                                ssb = p$ssb, rec = p$rec) }))
  res[, scenario := factor(scenario, levels = names(sc))]
  ref <- if(type == "cumulative") res[scenario == utils::tail(names(sc), 1)] else res[scenario == "all fleets"]
  imp <- if(type == "cumulative") rbindlist(lapply(seq_along(grp), function(k)
    data.table(year = yrs, group = names(grp)[k],
               impact = res[scenario == names(sc)[k]]$ssb - res[scenario == names(sc)[k + 1]]$ssb)))
  else rbindlist(lapply(names(grp), function(g)
    data.table(year = yrs, group = g, impact = res[scenario == paste0("-", g)]$ssb - ref$ssb)))
  b0d <- res[scenario == "no fishing"]$ssb
  imp[, prop := impact/(b0d[match(year, yrs)] - ref$ssb[match(year, yrs)])]
  imp[, group := factor(group, levels = names(grp))]
  chk <- data.table(year = yrs, ss3 = ssb3, projected = ref$ssb)
  dz <- as.data.table(ss3rep$Dynamic_Bzero)
  if(nrow(dz) && "SSB_nofishing" %in% names(dz))
    chk[, ss3_dynB0 := dz[match(yrs, Yr)]$SSB_nofishing][, dynB0 := b0d]
  flqs <- long2flqs(res[, .(name = as.character(scenario), age = "all", year, unit = "unique", season = "all",
                            area = "unique", iter = "1", data = ssb)], units = ss3rep$SpawnOutputUnits %||% "NA",
                    order = names(sc))
  structure(list(ssb = res, impact = imp, flqs = flqs, check = chk, groups = grp, type = type, init = init),
            class = "ss3impact")
}

#' @export
print.ss3impact <- function(x, ...){
  cat("Partial impact (Wang et al. 2009),", x$type, "- groups:", paste(names(x$groups), collapse = ", "), "\n")
  e <- x$check[, max(abs(projected/ss3 - 1))]
  cat(sprintf("max relative difference SS3 vs projected SSB (all fleets): %.2g\n", e))
  if("ss3_dynB0" %in% names(x$check))
    cat(sprintf("max relative difference SS3 vs projected dynamic B0 (init = \"%s\"): %.2g\n", x$init,
                x$check[, max(abs(dynB0/ss3_dynB0 - 1), na.rm = TRUE)]))
  ly <- max(x$impact$year)
  cat("impact in", ly, ":\n"); print(x$impact[year == ly, .(group, impact, prop)])
  invisible(x)
}

#' Plot the partial impact of fleets
#'
#' Plots \code{\link{ss3_impact}()} as in Wang et al. (2009): the SSB (grey) with the SSB lost to each fleet group
#' stacked on top, up to the dynamic unfished SSB (dashed line), or the share of the impact by group, or the SSB
#' trajectories of the scenarios.
#'
#' @param x object from \code{\link{ss3_impact}()}
#' @param type \code{"area"} (default), \code{"prop"} (share of the total impact) or \code{"line"} (SSB by scenario)
#' @param ss3 add the SS3 dynamic B0 (dotted), when reported
#' @return ggplot
#' @export
ss3plot_impact <- function(x, type = c("area", "prop", "line"), ss3 = TRUE){
  type <- match.arg(type)
  grp <- levels(x$impact$group)
  cols <- stats::setNames(ss3cols(length(grp)), grp)
  if(type == "line"){
    d <- copy(x$ssb)
    lc <- stats::setNames(c("black", ss3cols(nlevels(d$scenario) - 1)), levels(d$scenario))
    return(ggplot2::ggplot(d, ggplot2::aes(year, ssb, colour = scenario)) + ggplot2::geom_line(linewidth = 0.7) +
             ggplot2::scale_colour_manual(values = lc) + ggplot2::expand_limits(y = 0) +
             ggplot2::labs(x = "Year", y = "SSB", colour = NULL) + scale_x_year() + theme_ss3())
  }
  if(x$type != "cumulative" && type == "area")
    warning("impacts of type \"remove\" do not add up; consider type = \"line\" or \"prop\"", call. = FALSE)
  if(type == "prop")
    return(ggplot2::ggplot(x$impact, ggplot2::aes(year, prop, fill = group)) +
             ggplot2::geom_area(position = "fill", colour = "white", linewidth = 0.1) +
             ggplot2::scale_fill_manual(values = cols) +
             ggplot2::labs(x = "Year", y = "Share of the impact on SSB", fill = NULL) + scale_x_year() + theme_ss3())
  ref <- if(x$type == "cumulative") x$ssb[scenario == utils::tail(levels(scenario), 1)] else x$ssb[scenario == "all fleets"]
  d <- rbind(data.table(year = ref$year, group = "SSB", v = ref$ssb), x$impact[, .(year, group, v = impact)])
  d[, group := factor(group, levels = rev(c("SSB", grp)))]
  p <- ggplot2::ggplot(d, ggplot2::aes(year, v, fill = group)) +
    ggplot2::geom_area(position = "stack", colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = c(SSB = "grey65", cols), breaks = c("SSB", grp)) +
    ggplot2::geom_line(data = x$ssb[scenario == "no fishing"], ggplot2::aes(year, ssb), inherit.aes = FALSE,
                       linetype = 2, linewidth = 0.6)
  if(ss3 && "ss3_dynB0" %in% names(x$check))
    p <- p + ggplot2::geom_line(data = x$check, ggplot2::aes(year, ss3_dynB0), inherit.aes = FALSE, linetype = 3,
                                linewidth = 0.6, colour = "grey20")
  p + ggplot2::labs(x = "Year", y = "SSB", fill = NULL,
                    caption = paste0("dashed: SSB without fishing (dynamic B0)",
                                     if(ss3 && "ss3_dynB0" %in% names(x$check)) "; dotted: SS3 dynamic B0")) +
    scale_x_year() + theme_ss3()
}
