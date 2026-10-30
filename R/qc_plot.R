# Four views of a chain result, base graphics only. Every function draws and
# hands back the numbers it drew, invisibly, so the same table can go into a
# report. One colour per level, the same in every plot.

# the colour of each level code; interpolation codes are grey
qc_level_colours <- function() {
  c("1" = "#d62728", "2" = "#ff7f0e", "3" = "#9467bd", "4" = "#1f77b4", "5" = "#8c564b",
    "6" = "#2ca02c", "7" = "#e377c2", "8" = "#17becf", "9" = "#bcbd22", "50" = "grey70")
}

# level code per cell: 0 unobjected, N objected by level N (a refilled removal
# 50 + N counts as N), NA where nothing was observed or a plain gap was filled
qc_level_codes <- function(res) {
  F <- zoo::coredata(res$qc_data_flagged)
  L <- ifelse(!is.na(F) & F >= 50, F - 50, F)
  L[!is.na(L) & L == 0 & !is.na(F) & F == 50] <- NA
  L
}

# level labels of the chain this result belongs to, keyed by code
qc_level_labels <- function(res) {
  chain <- if (!is.null(res$qc_info$dataset_humidity)) "rh" else "t"
  keys  <- grep(paste0("^", chain, "[0-9]_"), names(res$qc_info), value = TRUE)
  codes <- as.integer(sub(paste0("^", chain, "([0-9])_.*$"), "\\1", keys))
  stats::setNames(keys, codes)
}

# the result must be a chain result
qc_check_result <- function(res) {
  if (!is.list(res) || !all(c("qc_data", "qc_data_flagged", "qc_info") %in% names(res)))
    stop("res must be a chain result list with qc_data, qc_data_flagged and qc_info.")
}

#' Flags per day, one line per level
#'
#' When did the chain object, and which level did it? A front shows as a
#' spike across the network, a failed sensor as a step that stays until the
#' end of the record.
#'
#' @param res A chain result list.
#' @param station One station name to plot alone; NULL for the whole network.
#' @param main Plot title.
#'
#' @return Invisibly, a data frame with day, level and n (flags that day).
#'
#' @examples
#' \dontrun{
#' qc_plot_timeline(t_res)
#' qc_plot_timeline(t_res, station = "Log_B4")
#' }
#'
#' @import zoo
#' @export
qc_plot_timeline <- function(res, station = NULL, main = NULL) {
  qc_check_result(res)
  L <- qc_level_codes(res)
  if (!is.null(station)) {
    if (!station %in% colnames(L)) stop(sprintf("Station '%s' is not in the result.", station))
    L <- L[, station, drop = FALSE]
  }
  day <- as.Date(format(zoo::index(res$qc_data_flagged), "%Y-%m-%d"))
  days <- sort(unique(day))
  labels <- qc_level_labels(res)
  codes <- sort(unique(L[!is.na(L) & L > 0]))
  # flags per day and level; an empty matrix when nothing was objected to
  tab <- matrix(0, nrow = length(days), ncol = length(codes), dimnames = list(as.character(days), codes))
  for (j in seq_along(codes)) tab[, j] <- rowsum(as.numeric(rowSums(L == codes[j], na.rm = TRUE)), day)[, 1]
  out <- data.frame(day = rep(days, times = length(codes)), level = rep(labels[as.character(codes)], each = length(days)),
                    n = as.vector(tab), stringsAsFactors = FALSE)
  cols <- qc_level_colours()
  plot(days, if (length(codes)) tab[, 1] else rep(0, length(days)), type = "n",
       ylim = c(0, max(1, tab)), xlab = "", ylab = "flags per day",
       main = if (is.null(main)) paste("Flags per day", if (!is.null(station)) station else "") else main)
  grid(nx = NA, ny = NULL, col = "grey90")
  for (j in seq_along(codes))
    lines(days, tab[, j], col = cols[as.character(codes[j])], lwd = 2)
  if (length(codes))
    legend("topleft", legend = labels[as.character(codes)], col = cols[as.character(codes)],
           lwd = 2, bty = "n", cex = 0.8)
  invisible(out)
}

#' Flags per level, or per station
#'
#' `by = "level"`: one bar per level, with the share of the judged values it
#' objected to written above - so a small count at a level that saw only one
#' station reads differently from a small count at a level that saw them all.
#' `by = "station"`: one horizontal bar per station, stacked by level, sorted
#' by the share of its observed values that were objected to - the list of
#' conspicuous stations.
#'
#' @param res A chain result list.
#' @param by `"level"` or `"station"`.
#' @param main Plot title.
#'
#' @return Invisibly: for `by = "level"` a data frame with level, flagged,
#'   judged and share; for `by = "station"` the matrix of shares (levels in
#'   rows, stations in columns, in percent of the station's observed values).
#'
#' @examples
#' \dontrun{
#' qc_plot_levels(t_res)
#' qc_plot_levels(t_res, by = "station")
#' }
#'
#' @import zoo
#' @export
qc_plot_levels <- function(res, by = c("level", "station"), main = NULL) {
  qc_check_result(res)
  by <- match.arg(by)
  L <- qc_level_codes(res)
  labels <- qc_level_labels(res)
  cols <- qc_level_colours()
  if (by == "level") {
    s <- qc_summary(res)
    agg <- do.call(rbind, lapply(split(s, factor(s$level, levels = unique(s$level))), function(z)
      data.frame(level = z$level[1], flagged = sum(z$flagged), judged = sum(z$judged), stringsAsFactors = FALSE)))
    agg$share <- ifelse(agg$judged > 0, 100 * agg$flagged / agg$judged, NA)
    rownames(agg) <- NULL
    codes <- names(labels)[match(agg$level, labels)]
    op <- graphics::par(mar = c(11, 4, 3, 1)); on.exit(graphics::par(op))
    bp <- graphics::barplot(agg$flagged, names.arg = agg$level, las = 2, col = cols[codes],
                            ylab = "flagged values", cex.names = 0.75,
                            ylim = c(0, max(1, agg$flagged) * 1.15),
                            main = if (is.null(main)) "Flags per level" else main)
    graphics::text(bp, agg$flagged, pos = 3, cex = 0.75,
                   labels = ifelse(is.na(agg$share), "not judged", sprintf("%d (%.2f%% of judged)", agg$flagged, agg$share)))
    return(invisible(agg))
  }
  # by station: share of each station's observed values objected to, per level
  n_obs <- colSums(!is.na(L))
  codes <- sort(unique(L[!is.na(L) & L > 0]))
  M <- sapply(colnames(L), function(st) sapply(codes, function(k) 100 * sum(L[, st] == k, na.rm = TRUE) / max(1, n_obs[st])))
  M <- matrix(M, nrow = length(codes), dimnames = list(labels[as.character(codes)], colnames(L)))
  if (!length(codes)) {
    graphics::plot.new(); graphics::title(main = if (is.null(main)) "Flags per station" else main)
    graphics::text(0.5, 0.5, "no flags"); return(invisible(M))
  }
  ord <- order(colSums(M))
  op <- graphics::par(mar = c(4, 8, 3, 1)); on.exit(graphics::par(op))
  graphics::barplot(M[, ord, drop = FALSE], horiz = TRUE, las = 1, col = cols[as.character(codes)],
                    xlab = "% of the station's observed values", cex.names = 0.7,
                    main = if (is.null(main)) "Flags per station" else main)
  if (length(codes))
    graphics::legend("bottomright", legend = labels[as.character(codes)], fill = cols[as.character(codes)], bty = "n", cex = 0.7)
  invisible(M)
}

#' Station-by-day calendar of the flags
#'
#' One row per station, one column per day, the colour is the share of the
#' day's observed values that were objected to; days without observations
#' are grey. Network-wide weather days appear as vertical stripes, a broken
#' sensor as a horizontal bar, data gaps as grey blocks.
#'
#' @param res A chain result list.
#' @param level One level key (e.g. `"t6_spatial_consistency"`) to show that
#'   level's flags alone; NULL for all levels together.
#' @param main Plot title.
#'
#' @return Invisibly, the matrix of shares (days in rows, stations in
#'   columns, in percent; NA where nothing was observed).
#'
#' @examples
#' \dontrun{
#' qc_plot_calendar(t_res)
#' qc_plot_calendar(t_res, level = "t6_spatial_consistency")
#' }
#'
#' @import zoo
#' @export
qc_plot_calendar <- function(res, level = NULL, main = NULL) {
  qc_check_result(res)
  L <- qc_level_codes(res)
  labels <- qc_level_labels(res)
  if (!is.null(level)) {
    if (!level %in% labels) stop(sprintf("level must be one of: %s", paste(labels, collapse = ", ")))
    code <- as.integer(names(labels)[labels == level])
    hit <- L == code
  } else {
    hit <- L > 0
  }
  day <- as.Date(format(zoo::index(res$qc_data_flagged), "%Y-%m-%d"))
  days <- sort(unique(day))
  obs <- rowsum((!is.na(L)) * 1, day)
  nh  <- rowsum(ifelse(is.na(hit), 0, hit) * 1, day)
  M <- 100 * nh / obs
  M[obs == 0] <- NA
  M <- matrix(M, nrow = length(days), dimnames = list(as.character(days), colnames(L)))
  ns <- ncol(M)
  op <- graphics::par(mar = c(6, 8, 3, 1)); on.exit(graphics::par(op))
  # classes rather than a linear ramp, so that a single objected value in a
  # day is still visible and a fully objected day stands out
  brk <- c(-1, 0, 1, 5, 10, 25, 50, 100)
  pal <- c("white", "#fff1e0", "#fdd49e", "#fdbb84", "#fc8d59", "#e34a33", "#7f0000")
  graphics::image(as.numeric(days), seq_len(ns), M, col = pal, breaks = brk,
                  xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = if (is.null(main)) paste("Flagged share per day", if (!is.null(level)) level else "") else main)
  # grey where nothing was observed
  na <- which(is.na(M), arr.ind = TRUE)
  if (nrow(na)) graphics::rect(as.numeric(days)[na[, 1]] - 0.5, na[, 2] - 0.5,
                               as.numeric(days)[na[, 1]] + 0.5, na[, 2] + 0.5, col = "grey85", border = NA)
  graphics::axis(2, at = seq_len(ns), labels = colnames(M), las = 1, cex.axis = 0.6)
  at <- pretty(days)
  graphics::axis(1, at = as.numeric(at), labels = format(at, "%d.%m."), cex.axis = 0.8)
  # the colour key, below the plot
  graphics::legend("bottom", inset = c(0, -0.28), xpd = TRUE, bty = "n", cex = 0.7, horiz = TRUE,
                   legend = c("0", "< 1 %", "1-5", "5-10", "10-25", "25-50", "> 50 %", "no data"),
                   fill = c(pal, "grey85"))
  invisible(M)
}

#' One station's series with its flags
#'
#' The cleaned series as a line, the raw series behind it in grey when it is
#' given, and the objected cells as coloured points by level. Without the
#' raw series the objected cells are shown as ticks at the bottom, because
#' their values are no longer in the cleaned series.
#'
#' @param res A chain result list.
#' @param station Station name.
#' @param raw The series as it entered the chain (an xts with this station as
#'   a column), or NULL.
#' @param from,to Time window as anything `xts` understands (`"2025-07"`,
#'   `"2025-07-01/2025-07-15"` goes into `from`); NULL for everything.
#' @param main Plot title.
#'
#' @return Invisibly, a data frame with time, value (raw where given, else
#'   NA) and level for every objected cell in the window.
#'
#' @examples
#' \dontrun{
#' qc_plot_station(t_res, "Log_B4", raw = x)
#' qc_plot_station(t_res, "Log_B4", raw = x, from = "2025-08-05/2025-08-07")
#' }
#'
#' @import xts
#' @import zoo
#' @export
qc_plot_station <- function(res, station, raw = NULL, from = NULL, to = NULL, main = NULL) {
  qc_check_result(res)
  if (!station %in% colnames(res$qc_data)) stop(sprintf("Station '%s' is not in the result.", station))
  cln <- res$qc_data[, station]
  flg <- res$qc_data_flagged[, station]
  rw  <- if (!is.null(raw) && station %in% colnames(raw)) raw[, station] else NULL
  win <- if (!is.null(from)) paste0(from, if (!is.null(to)) paste0("/", to) else "") else NULL
  if (!is.null(win)) { cln <- cln[win]; flg <- flg[win]; if (!is.null(rw)) rw <- rw[win] }
  if (nrow(cln) == 0) stop("No data in that window.")
  lv <- as.numeric(zoo::coredata(flg)); lv <- ifelse(!is.na(lv) & lv >= 50, lv - 50, lv)
  hit <- !is.na(lv) & lv > 0
  labels <- qc_level_labels(res); cols <- qc_level_colours()
  ti <- zoo::index(cln)
  yr <- range(c(zoo::coredata(cln), if (!is.null(rw)) zoo::coredata(rw)), na.rm = TRUE)
  plot(ti, as.numeric(zoo::coredata(cln)), type = "n", ylim = yr, xlab = "", ylab = "",
       main = if (is.null(main)) station else main)
  grid(nx = NA, ny = NULL, col = "grey90")
  if (!is.null(rw)) graphics::lines(zoo::index(rw), as.numeric(zoo::coredata(rw)), col = "grey70")
  graphics::lines(ti, as.numeric(zoo::coredata(cln)), col = "black")
  vals <- if (!is.null(rw)) as.numeric(zoo::coredata(rw))[match(ti, zoo::index(rw))] else rep(NA_real_, length(ti))
  codes <- sort(unique(lv[hit]))
  for (k in codes) {
    sel <- hit & lv == k
    if (!is.null(rw)) graphics::points(ti[sel], vals[sel], col = cols[as.character(k)], pch = 16, cex = 0.6)
    else graphics::rug(as.numeric(ti[sel]), col = cols[as.character(k)], lwd = 2)
  }
  if (length(codes))
    graphics::legend("topleft", legend = labels[as.character(codes)], col = cols[as.character(codes)],
                     pch = if (!is.null(rw)) 16 else NA, lwd = if (!is.null(rw)) NA else 2, bty = "n", cex = 0.75)
  invisible(data.frame(time = ti[hit], value = vals[hit], level = unname(labels[as.character(lv[hit])]),
                       stringsAsFactors = FALSE))
}
