#' Convert a wide data frame to the xts object the QC chain expects
#' The input must have exactly one time column and one numeric column
#' per logger.
#'
#' @param x A data frame (or tibble): one time column plus one numeric column
#'   per logger. Logger columns need unique names; any number is allowed.
#' @param time_col The time column, required: its exact name (case-sensitive)
#'   or its column number.
#' @param time_format The format of the time stamps, as in `strptime()`.
#'   Ignored if the time column is already POSIXct.
#' @param tz Time zone of the time stamps.
#' @param verbose Print the conversion report.
#'
#' @return An `xts` object with a numeric matrix, one column per logger. It
#'   carries the conversion report as the xts attribute `qc_conversion`; the
#'   first QC level moves it into `qc_info$conversion_temperature` (or
#'   `_humidity`), so the finished result documents its own input.
#'
#' @examples
#' \dontrun{
#' raw <- read.csv("Biel_summer_2025_raw.csv", sep = ";", dec = ".",
#'                 fileEncoding = "UTF-8-BOM")
#' stopifnot(nrow(raw) == 13308)          # verify the read BEFORE converting
#' x <- qc_as_xts(raw, time_col = "time", time_format = "%d.%m.%Y %H:%M")
#' x <- qc_as_xts(raw, time_col = 1,      time_format = "%d.%m.%Y %H:%M")
#' }
#'
#' @import xts
#' @import zoo
#' @export
qc_as_xts <- function(x,
                      time_col,
                      time_format = "%Y-%m-%d %H:%M:%S",
                      tz = "UTC",
                      verbose = TRUE) {

  #-------------------------------------------------------------------------------
  # the table itself

  # only data frames (or tibbles) are accepted
  if (!is.data.frame(x))
    stop("x must be a data frame. Read the file in your script (e.g. read.csv), check it, then pass the result.")
  # There must be one time column plus at least one logger column
  if (ncol(x) < 2)
    stop("x needs a time column and at least one logger column.")
  # every column needs a unique, non-empty name
  nm <- names(x)
  if (anyNA(nm) || any(!nzchar(nm)))
    stop("Every column needs a non-empty name.")
  if (anyDuplicated(nm))
    stop(sprintf("Every column needs a unique name. Duplicated: %s",
                 paste(unique(nm[duplicated(nm)]), collapse = ", ")))

  #-------------------------------------------------------------------------------
  # the time column: named exactly by the caller, never searched for

  if (missing(time_col))
    stop("time_col is required: give the exact name or the number of the time column.")
  if (length(time_col) != 1)
    stop("time_col must be ONE column (name or a positive natural number).")
  if (is.numeric(time_col)) {
    # check whether the input is a positive natural number and within the range
    if (time_col != round(time_col) || time_col < 1 || time_col > ncol(x))
      stop(sprintf("time_col = %s is not a column number between 1 and %d.", time_col, ncol(x)))
    # assign the exact column name of the time column
    tc <- nm[time_col]
  } else if (is.character(time_col)) {
    # a column name: exact match, case-sensitive
    if (!time_col %in% nm)
      stop(sprintf("time_col '%s' not found. Columns are: %s", time_col, paste(nm, collapse = ", ")))
    # assign the exact column name of the time column
    tc <- time_col
  } else {
    stop("time_col must be a column name or a column number.")
  }

  #-------------------------------------------------------------------------------
  # the time stamps: parsed with exactly the format the caller gives

  # extract the time column
  tv <- x[[tc]]
  # true/false
  was_posixct <- inherits(tv, "POSIXct")
  if (!was_posixct)
    tv <- as.POSIXct(as.character(tv), format = time_format, tz = tz)
  # every stamp must exist - a POSIXct column can carry NA too, so this check
  # sits OUTSIDE the parsing branch; name the first one that is missing
  if (anyNA(tv)) {
    bad <- which(is.na(tv))
    stop(sprintf("%d time stamp(s) missing or not matching time_format '%s' (first: row %d, '%s').",
                 length(bad), time_format, bad[1], as.character(x[[tc]][bad[1]])))
  }
  # duplicated stamps are always an error.
  if (anyDuplicated(tv)) {
    d <- which(duplicated(tv))
    stop(sprintf("%d duplicated time stamp(s) (first: row %d, %s). If the raw stamps differ, time_format drops part of them.",
                 length(d), d[1], format(tv[d[1]])))
  }

  #-------------------------------------------------------------------------------
  # the logger columns: everything else, and all of it numeric

  # take all column names except the time column
  loggers <- setdiff(nm, tc)
  # check whether they are numeric
  ok <- vapply(x[loggers], function(v) is.numeric(v) || (is.logical(v) && all(is.na(v))), logical(1))
  if (!all(ok))
    stop(sprintf("Logger column(s) not numeric: %s. Decimal comma? Use read.csv(dec = \",\"). A text column? Remove it in your script.",
                 paste(loggers[!ok], collapse = ", ")))
  # plain numeric matrix, one column per logger
  m <- as.matrix(x[loggers])
  storage.mode(m) <- "double"

  #-------------------------------------------------------------------------------
  # report and hand the xts on to the first QC level

  # xts() sorts by time; say so if the input was not in order (e.g. appended files)
  resorted <- is.unsorted(tv)
  # the xts the chain works on
  out <- xts::xts(m, order.by = tv)
  # loggers without a single value (installed later, dead, wrong column)
  empty <- colnames(m)[colSums(!is.na(m)) == 0]

  # the standard record, in the same shape as every QC level's qc_info entry
  rec <- list(time_col      = tc,
              time_format   = if (was_posixct) "POSIXct (taken as is)" else time_format,
              tz            = tz,
              n_time        = nrow(out),
              n_loggers     = ncol(out),
              loggers       = colnames(out),
              time_start    = zoo::index(out)[1],
              time_end      = zoo::index(out)[nrow(out)],
              pct_missing   = round(100 * mean(is.na(m)), 1),
              empty_loggers = empty,
              resorted      = resorted)
  # the record travels WITH the xts: the first QC level moves it into
  # qc_info$conversion_<what>, next to qc_info$dataset_<what>
  xts::xtsAttributes(out) <- list(qc_conversion = rec)

  # the report, so the result of the read can be checked before the QC starts
  if (isTRUE(verbose)) {
    message(sprintf("qc_as_xts: %d time steps x %d loggers, %s to %s",
                    rec$n_time, rec$n_loggers, format(rec$time_start), format(rec$time_end)))
    message(sprintf("  time column '%s', format %s, tz %s", tc,
                    if (was_posixct) "POSIXct (taken as is)" else sprintf("'%s'", time_format), tz))
    message(sprintf("  missing values: %.1f%%", rec$pct_missing))
    if (length(empty))
      message("  empty loggers (not a single value): ", paste(empty, collapse = ", "))
    if (resorted)
      message("  rows were not in time order - sorted by time")
  }
  out
}
