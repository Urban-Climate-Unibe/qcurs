#' Temperature QC Level 5: climatic outliers of the daily extremes
#'
#' A station-specific, month-specific plausibility bound from the station's
#' own history. Per station and calendar month, all years pooled, the daily
#' maxima and the daily minima are collected; a day whose maximum exceeds
#' Q3 + m * IQR of the pooled daily maxima, or whose minimum falls below
#' Q1 - m * IQR of the pooled daily minima, is objected to - on the cells
#' beyond that bound, not the whole day. Only these outer tails are tested:
#' a cool summer day or a mild winter night is weather, not a sensor fault.
#'
#' The published chain applies the IQR rule to the 10-minute values, where
#' the bounds come out 25 K beyond the quartiles and never bite; on the
#' daily extremes with the published m = 4 they still sit 25 K out. With
#' m = 2.5 the bounds lie about 15 K above Q3 and 11 K below Q1 - Bern in June
#' 41 instead of the 46 of level 2, January -16 instead of -24 - and on six
#' years of three MeteoSwiss stations objected to four values in all: two
#' December nights at Thun, 0.9 K below its bound.
#'
#' Two guards apply: a station-month needs `min_days` days with at least
#' `min_obs_per_day` values, and - by default - its sample must span more
#' than one calendar year. A single-season campaign is therefore refused,
#' station-month by station-month, and says so. Set
#' `require_multi_year = FALSE` to override.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param ext_lim_factor Multiplier m on the IQR of the daily extremes.
#' @param min_days Minimum number of judged days a station-month needs.
#' @param min_obs_per_day Valid values a day needs to contribute an extreme
#'   (72 = half of a 10-minute day).
#' @param require_multi_year Only judge a station-month whose sample spans
#'   more than one calendar year.
#' @param verbose Report the tally and the refusals.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info. The
#'   refused station-months are listed in
#'   `qc_info$t5_climatic_outliers$refused`, the bounds that were applied in
#'   `qc_info$t5_climatic_outliers$bounds` (one row per judged
#'   station-month).
#'
#' @examples
#' \dontrun{
#' res <- T_QC_5_climatic_outliers(res)                             # honest self-skip
#' res <- T_QC_5_climatic_outliers(res, require_multi_year = FALSE) # force it anyway
#' }
#'
#' @import xts
#' @import zoo
#' @export
T_QC_5_climatic_outliers <- function(input,
                                     ext_lim_factor = 2.5,
                                     min_days = 60,
                                     min_obs_per_day = 72,
                                     require_multi_year = TRUE,
                                     verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "temperature", level = "t5_climatic_outliers")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # validate the parameters, because the caller may set them freely

  # the IQR multiplier must be positive
  if (!is.numeric(ext_lim_factor) || length(ext_lim_factor) != 1 || ext_lim_factor <= 0)
    stop("ext_lim_factor must be a positive number.")
  # quartiles from fewer than ~20 days are noise, not a climatology
  if (!is.numeric(min_days) || length(min_days) != 1 || min_days < 20)
    stop("min_days must be at least 20.")
  # an extreme from a handful of values is not a daily extreme
  if (!is.numeric(min_obs_per_day) || length(min_obs_per_day) != 1 || min_obs_per_day < 2)
    stop("min_obs_per_day must be at least 2.")

  #-------------------------------------------------------------------------------
  # Perform QC Level 5

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # calendar day of every time step, in the time zone of the index
  day <- format(zoo::index(x), "%Y-%m-%d")
  # the row numbers of every day, once, in the order of the record
  day_rows <- split(seq_along(day), factor(day, levels = unique(day)))
  # calendar month and year of every day, for the pooling and the guards
  day_mon <- as.numeric(substr(names(day_rows), 6, 7))
  day_yr  <- as.numeric(substr(names(day_rows), 1, 4))
  # count across all stations
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # per-station coverage: values on days of judged station-months
  n_judged <- stats::setNames(integer(ncol(X)), colnames(X))
  # which station-months were refused, so a zero-flag run is not mistaken for a clean one
  refused <- character(0)
  # the bounds that were applied, one row per judged station-month
  bounds <- list()

  # iterate over all stations (columns): a climatology is site-specific
  for (s in colnames(X)) {
    # extract data vector of this station
    v <- X[, s]
    # daily maxima and minima; NA where the day is too sparse
    ext <- vapply(day_rows, function(rows) {
      z <- v[rows]
      if (sum(!is.na(z)) < min_obs_per_day) return(c(NA_real_, NA_real_))
      range(z, na.rm = TRUE)
    }, numeric(2))
    tmin <- ext[1, ]; tmax <- ext[2, ]
    # every calendar month present in the record
    for (mm in sort(unique(day_mon))) {
      # this station's judged days in this month, all years pooled
      sel <- day_mon == mm & !is.na(tmax)
      # sample size of the pooled month, in days
      nd <- sum(sel)
      # over how many distinct years the sample spreads
      years_span <- length(unique(day_yr[sel]))
      # guards: enough days AND (if required) a genuine multi-year climatology
      if (nd < min_days || (require_multi_year && years_span < 2)) {
        # remember the refusal instead of silently doing nothing
        if (nd > 0) refused <- c(refused, sprintf("%s/%02d(days=%d,y=%d)", s, mm, nd, years_span))
        next
      }
      # the outer tails: hot days by the maxima, cold nights by the minima
      q_hi <- stats::quantile(tmax[sel], c(0.25, 0.75))
      hi <- q_hi[2] + ext_lim_factor * (q_hi[2] - q_hi[1])
      q_lo <- stats::quantile(tmin[sel], c(0.25, 0.75))
      lo <- q_lo[1] - ext_lim_factor * (q_lo[2] - q_lo[1])
      bounds[[length(bounds) + 1]] <- data.frame(station = s, month = mm, days = nd,
                                                 tmax_upper = unname(hi), tmin_lower = unname(lo),
                                                 stringsAsFactors = FALSE)
      # this station-month is judged, whatever the verdict
      rows_m <- unlist(day_rows[sel], use.names = FALSE)
      n_judged[s] <- n_judged[s] + sum(!is.na(v[rows_m]))
      # the offending days, and on them the cells beyond the bound
      hit <- rep(FALSE, length(v))
      for (d in which(sel & (tmax > hi | tmin < lo))) {
        rows <- day_rows[[d]]
        hit[rows] <- !is.na(v[rows]) & ((tmax[d] > hi & v[rows] > hi) | (tmin[d] < lo & v[rows] < lo))
      }

      # combine the verdict with THIS station's column only
      mask <- hit & (is.na(previous_flag[, s]) | previous_flag[, s] == 0)
      # how many cells this level objects to in this station-month
      n_found <- sum(mask)

      # apply only if something was found
      if (n_found > 0) {
        # blank the outliers so later levels never see them
        X[mask, s] <- NA
        # record this level's code (5 = level 5, fixed by convention)
        previous_flag[mask, s] <- 5
        # add the number of new flags to the counters
        n_station[s] <- n_station[s] + n_found
        n_total <- n_total + n_found
      }
    }
  }

  # write the matrices back into the xts shells, keeping index and column names
  if (n_total > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }
  # the bounds as one table, empty when nothing was judged
  bounds <- if (length(bounds)) do.call(rbind, bounds) else
    data.frame(station = character(0), month = numeric(0), days = integer(0),
               tmax_upper = numeric(0), tmin_lower = numeric(0), stringsAsFactors = FALSE)

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # report so a zero-hit run is visibly a run, not a skip - and name the refusals
  if (isTRUE(verbose)) {
    message(sprintf("T5 climatic outliers (daily extremes, m=%g): %d flagged, %d station-months judged, %d refused",
                    ext_lim_factor, n_total, nrow(bounds), length(refused)))
    if (length(refused) > 0 && length(refused) <= 12)
      message("  refused: ", paste(refused, collapse = ", "))
    else if (length(refused) > 12)
      message("  refused (first 12): ", paste(utils::head(refused, 12), collapse = ", "), ", ...")
  }

  # write the updated matrices back and append this level under its own name
  input$qc_data                      <- x
  input$qc_data_flagged              <- flg
  input$qc_info$t5_climatic_outliers <- list(n_flagged            = n_total,
                                             n_flagged_by_station = n_station,
                                             n_judged_by_station  = n_judged,
                                             ext_lim_factor       = ext_lim_factor,
                                             min_days             = min_days,
                                             min_obs_per_day      = min_obs_per_day,
                                             require_multi_year   = require_multi_year,
                                             refused              = refused,
                                             bounds               = bounds)
  input
}
