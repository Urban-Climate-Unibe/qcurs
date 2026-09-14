#' Convert a data frame to the xts object the QC chain expects
#'
#' Handles the one conversion that is genuinely ambiguous and therefore worth
#' centralising: a wide data frame with one time column and one column per
#' station. Building the xts by hand with `xts(x = df, order.by = df$time)`
#' keeps the time column inside the matrix and yields a CHARACTER matrix, in
#' which every comparison is lexical - the single most common way to make a QC
#' run look fine and be wrong.
#'
#' The time column is found by name if you give one, otherwise by matching the
#' column names against a list of common spellings, case- and
#' punctuation-insensitively ("time", "TIMESTAMP", "Date.Time..GMT.02.00",
#' "Zeitstempel", ...). Ambiguity is never resolved by guessing: if several
#' columns qualify, or if none does, the function stops and names the
#' candidates so you can pass `time_col` explicitly. Loggers that export the
#' date and the clock time in two columns are supported by passing both, e.g.
#' `time_col = c("Date", "Time")`.
#'
#' Time stamps are parsed with the first candidate format that (a) parses
#' EVERY stamp and (b) does not collapse distinct strings onto identical
#' stamps. Condition (b) is not cosmetic: `as.POSIXct()` silently ignores the
#' trailing part of a string, so the date-only format "%Y-%m-%d" parses
#' "2025-06-01 00:10" without complaint and maps a whole day onto midnight.
#' Date-only formats are therefore tried last AND must survive the collapse
#' check. Numeric columns are accepted as Unix epoch seconds (or
#' milliseconds) if they fall in a plausible range.
#'
#' Both layouts are accepted. WIDE means one column per station, which is what
#' the QC chain works on internally - the spatial levels compare stations
#' column-wise at one time step, so the chain cannot operate on long data.
#' LONG ("tidy") means one row per observation with a station column and a
#' value column; it is the better storage and interchange format and is
#' pivoted to wide here. The pivot is strict: a repeated (time, station) pair
#' is an error, never a silent overwrite, and every missing combination
#' becomes NA so the time grid stays complete.
#'
#' File reading is deliberately NOT part of this function. Reading a CSV
#' involves separator, decimal mark, encoding and time format decisions whose
#' failures are silent (a wrong encoding can truncate the file mid-way and
#' still "succeed"), so it belongs in the caller's script where the result can
#' be inspected, not inside the QC chain.
#'
#' @param x A data frame (or matrix) with one time column and one numeric
#'   column per station, or an `xts` object, which is returned unchanged.
#' @param time_col Which column holds the time. Either character (column
#'   name, matched case- and punctuation-insensitively) or numeric (column
#'   POSITION, 1-based, for files whose header is unusable). Give two values
#'   for loggers that split date and clock time, e.g. `c("Date", "Time")` or
#'   `c(1, 2)`. NULL (default) auto-detects by name. Note that `c("Date", 2)`
#'   does NOT mix the two: R coerces it to character and the "2" is then
#'   looked up as a column NAME.
#' @param verbose Report the detected column, the winning format and the
#'   interpretation of a numeric time column.
#' @param station_col Long ("tidy") input only: the column holding the station
#'   ID. NULL (default) auto-detects among the usual names ("station", "id",
#'   "sensor", "logger", "site"); the wide layout is assumed when none is
#'   found. Character name or numeric position, like `time_col`.
#' @param value_col Long input only: the column holding the measured value.
#'   NULL auto-detects when exactly one column is left over.
#' @param tz Time zone used when parsing character time stamps.
#' @param formats Candidate time formats, tried in order. Date-time formats
#'   must precede date-only ones.
#'
#' @return An `xts` object with a numeric matrix and one column per station.
#'
#' @examples
#' \dontrun{
#' raw <- read.csv("Bern_summer_2025_raw.csv", sep = ";", dec = ".",
#'                 fileEncoding = "UTF-8-BOM")
#' stopifnot(nrow(raw) == 13308)          # verify the read BEFORE converting
#' x <- qc_as_xts(raw)                    # finds "time", "TIMESTAMP", ...
#' x <- qc_as_xts(raw, time_col = c("Date", "Time"))   # split date/time columns
#' x <- qc_as_xts(tidy)                                # long: time, station, temperature
#' x <- qc_as_xts(tidy, station_col = "station", value_col = "temperature")
#' res <- T_QC_1_gross_error(x)
#' }
#'
#' @import xts
#' @import zoo
#' @export
qc_as_xts <- function(x, time_col = NULL, station_col = NULL, value_col = NULL,
                      tz = "UTC",
                      formats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
                                  "%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S",
                                  "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M",
                                  "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
                                  "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                                  "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
                                  "%Y%m%d%H%M%S",
                                  "%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y"),
                      verbose = TRUE) {
  
  #-------------------------------------------------------------------------------
  # the easy cases first
  
  # an xts is already what we want: hand it back untouched
  if (inherits(x, "xts")) return(x)
  # a bare path is refused on purpose - see the note above
  if (is.character(x) && length(x) == 1)
    stop("qc_as_xts() does not read files. Read the file in your script, check the result, then pass the data frame.")
  # everything else must be rectangular
  if (!is.data.frame(x) && !is.matrix(x)) stop("x must be a data frame, a matrix or an xts object.")
  x <- as.data.frame(x)
  if (ncol(x) < 2) stop("Need at least two columns (a time column and a data column).")
  # duplicated column names would silently drop stations further down: setdiff() de-duplicates,
  # so a second column called "Log_1" disappears without a word. Refuse instead.
  if (anyDuplicated(names(x)))
    stop(sprintf("Duplicated column names: %s. Every station needs a unique column name.",
                 paste(unique(names(x)[duplicated(names(x))]), collapse = ", ")))
  
  #-------------------------------------------------------------------------------
  # locate the time column: by name if given, otherwise by matching common spellings
  
  # strip a UTF-8 BOM, lowercase, drop every separator: "Date.Time..GMT." -> "datetimegmt"
  normalise <- function(v) gsub("[^a-z0-9]", "", tolower(sub("^\ufeff", "", v)))
  nm  <- names(x)
  key <- normalise(nm)
  
  #-------------------------------------------------------------------------------
  # check which column is the time column when the input is NOT NULL
  
  if (!is.null(time_col)) {
    # at most two columns can define one time index (date + clock time)
    if (length(time_col) < 1 || length(time_col) > 2) stop("time_col must name one column, or two for a split date/time layout.")
    
    #-------------------------------------------------------------------------------
    # if the time_column-input is numeric (the caller gives the column position)
    
    if (is.numeric(time_col)) {
      # check whether the input is a natural number
      if (any(time_col != round(time_col))) stop("A numeric time_col must be a whole column number.")
      # check whether the number is within the range
      if (any(time_col < 1 | time_col > ncol(x))) stop(sprintf("time_col %s is out of range: the data have %d columns.",paste(time_col[time_col < 1 | time_col > ncol(x)], collapse = ", "), ncol(x)))
      # check whether the same column is used more than once
      if (anyDuplicated(time_col)) stop("time_col names the same column twice.")
      # this is the column which is used as time_col
      tc <- nm[time_col]
      if (isTRUE(verbose))
        message(sprintf("qc_as_xts: using column %s ('%s') as the time column.",
                        paste(time_col, collapse = "+"), paste(tc, collapse = "', '")))
      
      #-------------------------------------------------------------------------------
      # if the time_column-input is a character
      
    } else if (is.character(time_col)) {
      # character: match the name case- and punctuation-insensitively, but insist it exists
      hit <- match(normalise(time_col), key)
      if (anyNA(hit)) {
        missing_nm <- time_col[is.na(hit)]
        # a digits-only string is almost always a position passed the wrong way
        extra <- if (all(grepl("^[0-9]+$", missing_nm)))
          " (for a column POSITION pass a number, not a string: time_col = 1)" else ""
        stop(sprintf("time_col '%s' not found%s. Available columns: %s",
                     paste(missing_nm, collapse = "', '"), extra, paste(nm, collapse = ", ")))
      }
      if (anyDuplicated(hit)) stop("time_col names the same column twice.")
      # this is the column which is used as time_col
      tc <- nm[hit]
      
      #-------------------------------------------------------------------------------
      # if the time_column-input is neither numeric nor a character
      
    } else {
      stop("time_col must be a column name (character) or a column position (numeric).")
    }
    
    #-------------------------------------------------------------------------------
    # if no time_col was given at all, find it among the common spellings
    
  } else {
    # the spellings loggers and portals actually use
    known <- c("time", "times", "timestamp", "timestamps", "timestmp",
               "datetime", "date", "datum", "zeit", "zeitstempel", "uhrzeit",
               "utc", "utctime", "localtime", "obstime", "observationtime",
               "measurementtime", "recordtime", "referencetimestamp", "stamp", "dt", "ts")
    # exact matches first
    cand <- which(key %in% known)
    # then prefix matches, which catch "datetimegmt0200", "timestampiso", ...
    if (!length(cand))
      cand <- grep("^(timestamp|datetime|datum|zeit|time|date)", key)
    if (!length(cand)) stop(sprintf("No time column found. Pass time_col explicitly. Available columns: %s",paste(nm, collapse = ", ")))
    if (length(cand) > 1) {
      # a date column AND a clock column is a known logger layout - say how to fix it
      is_date <- key[cand] %in% c("date", "datum")
      is_time <- key[cand] %in% c("time", "zeit", "uhrzeit")
      if (sum(is_date) == 1 && sum(is_time) == 1)
        stop(sprintf("Found a separate date and time column ('%s', '%s'). Pass time_col = c('%s', '%s').",
                     nm[cand][is_date], nm[cand][is_time], nm[cand][is_date], nm[cand][is_time]))
      # otherwise refuse to guess which one is the real index
      stop(sprintf("Several columns could be the time column: %s. Pass time_col explicitly.",
                   paste(nm[cand], collapse = ", ")))
    }
    tc <- nm[cand]
    if (isTRUE(verbose)) message(sprintf("qc_as_xts: using '%s' as the time column.", tc))
  }
  
  #-------------------------------------------------------------------------------
  # decide whether this is a long (tidy) or a wide layout
  
  # small helper: accept a column name (character) or a position (numeric), like time_col
  resolve <- function(v, label) {
    if (is.numeric(v)) {
      if (length(v) != 1 || v != round(v) || v < 1 || v > ncol(x))
        stop(sprintf("%s must be one valid column position (1..%d).", label, ncol(x)))
      return(nm[v])
    }
    if (!is.character(v) || length(v) != 1) stop(sprintf("%s must be one column name or position.", label))
    hit <- match(normalise(v), key)
    if (is.na(hit)) stop(sprintf("%s '%s' not found. Available columns: %s", label, v, paste(nm, collapse = ", ")))
    nm[hit]
  }
  
  if (!is.null(station_col)) {
    # the caller says it is long
    sc <- resolve(station_col, "station_col")
  } else {
    # auto-detect: the names a station column actually carries
    st_known <- c("station", "stationid", "stationname", "id", "sensor", "sensorid",
                  "logger", "loggerid", "site", "siteid", "name")
    hit <- which(key %in% st_known & !(nm %in% tc))
    # several candidates would make the pivot ambiguous - refuse instead of guessing
    if (length(hit) > 1)
      stop(sprintf("Several columns could be the station column: %s. Pass station_col explicitly.",
                   paste(nm[hit], collapse = ", ")))
    sc <- if (length(hit) == 1) nm[hit] else NULL
  }
  
  if (!is.null(sc)) {
    # long layout: find the value column among what is left
    rest <- setdiff(nm, c(tc, sc))
    if (!is.null(value_col)) {
      vc <- resolve(value_col, "value_col")
      if (vc %in% c(tc, sc)) stop("value_col must differ from time_col and station_col.")
    } else if (length(rest) == 1) {
      vc <- rest
    } else {
      stop(sprintf("Long layout detected (station column '%s') but %d possible value columns: %s. Pass value_col explicitly.",
                   sc, length(rest), paste(rest, collapse = ", ")))
    }
    if (isTRUE(verbose))
      message(sprintf("qc_as_xts: long layout - station '%s', value '%s'.", sc, vc))
  } else {
    vc <- NULL
  }
  
  #-------------------------------------------------------------------------------
  # parse the time stamps
  
  # one column, or two that get pasted together with a space
  tv <- if (length(tc) == 2) paste(as.character(x[[tc[1]]]), as.character(x[[tc[2]]])) else x[[tc]]
  
  if (!inherits(tv, "POSIXct")) {
    
    #-------------------------------------------------------------------------------
    # a numeric index: epoch seconds, epoch milliseconds, or a packed YYYYMMDDHHMMSS
    
    if (is.numeric(tv)) {
      rng <- range(tv, na.rm = TRUE)
      if (rng[1] > 9.4e8 && rng[2] < 4.1e9) {
        tv <- as.POSIXct(tv, origin = "1970-01-01", tz = tz)
        if (isTRUE(verbose)) message("qc_as_xts: time column read as Unix epoch seconds.")
      } else if (rng[1] > 9.4e11 && rng[2] < 4.1e12) {
        tv <- as.POSIXct(tv / 1000, origin = "1970-01-01", tz = tz)
        if (isTRUE(verbose)) message("qc_as_xts: time column read as Unix epoch milliseconds.")
      } else if (rng[1] > 1.9e13 && rng[2] < 2.2e13) {
        # read.csv turns a bare YYYYMMDDHHMMSS stamp into a number; recover it as text
        tv <- as.POSIXct(sprintf("%.0f", tv), format = "%Y%m%d%H%M%S", tz = tz)
        if (anyNA(tv)) stop(sprintf("Numeric time column '%s' looks like YYYYMMDDHHMMSS but does not parse.",
                                    paste(tc, collapse = "+")))
        if (isTRUE(verbose)) message("qc_as_xts: time column read as YYYYMMDDHHMMSS.")
      } else {
        stop(sprintf("Numeric time column '%s' is neither a plausible Unix epoch nor YYYYMMDDHHMMSS. Convert it in your script.",
                     paste(tc, collapse = "+")))
      }
      
      #-------------------------------------------------------------------------------
      # character stamps: try the candidate formats in order
      
    } else {
      s <- trimws(as.character(tv))
      # a missing stamp cannot become an xts index. Say that precisely, instead of
      # letting every format fail and blaming the format list.
      empty <- is.na(s) | !nzchar(s)
      if (any(empty))
        stop(sprintf("%d of %d time stamps in '%s' are missing or empty (first at row %d). Drop or repair those rows in your script.",
                     sum(empty), length(s), paste(tc, collapse = "+"), which(empty)[1]))
      # how many DISTINCT strings the column contains - the yardstick for the collapse check
      n_uniq <- length(unique(s))
      parsed <- NULL; used <- NA_character_
      for (f in formats) {
        cand <- as.POSIXct(s, format = f, tz = tz)
        # (a) the format must parse every single stamp ...
        if (anyNA(cand)) next
        # ... and (b) must not map distinct strings onto identical stamps.
        # Without this, "%Y-%m-%d" happily "parses" "2025-06-01 00:10" to midnight.
        if (length(unique(cand)) < n_uniq) next
        parsed <- cand; used <- f; break
      }
      if (is.null(parsed))
        stop(sprintf("Could not parse column '%s' as time with any of: %s",
                     paste(tc, collapse = "+"), paste(formats, collapse = ", ")))
      # a purely numeric date like 06/07/2025 fits day-first AND month-first. Both parse
      # everything, so the loop simply takes the earlier one - warn rather than let a
      # six-month shift pass unnoticed.
      if (grepl("%d[./]%m|%m[./]%d", used)) {
        twin <- if (grepl("%d([./])%m", used)) sub("%d([./])%m", "%m\\1%d", used) else sub("%m([./])%d", "%d\\1%m", used)
        alt <- as.POSIXct(s, format = twin, tz = tz)
        if (!anyNA(alt) && !isTRUE(all.equal(as.numeric(alt), as.numeric(parsed))))
          warning(sprintf("Ambiguous time stamps: '%s' and '%s' both parse this column. Used '%s' (%s, not %s). Pass `formats` explicitly if that is wrong.",
                          used, twin, used, format(parsed[1]), format(alt[1])))
      }
      if (isTRUE(verbose)) message(sprintf("qc_as_xts: time format '%s'.", used))
      tv <- parsed
    }
  }
  
  #-------------------------------------------------------------------------------
  # long layout: pivot one row per observation into one column per station
  
  if (!is.null(sc)) {
    st <- as.character(x[[sc]])
    if (anyNA(st) || any(!nzchar(trimws(st))))
      stop(sprintf("Station column '%s' has missing or empty entries (first at row %d).",
                   sc, which(is.na(st) | !nzchar(trimws(st)))[1]))
    st <- trimws(st)
    # the target grid: every time stamp once, every station once, both sorted
    ut <- sort(unique(tv))
    us <- sort(unique(st))
    ri <- match(tv, ut)
    ci <- match(st, us)
    # a repeated (time, station) pair would silently overwrite - that is a data error
    dup <- duplicated(cbind(ri, ci))
    if (any(dup))
      stop(sprintf("%d repeated (time, station) pairs, e.g. %s at %s. Aggregate or de-duplicate before converting.",
                   sum(dup), st[which(dup)[1]], format(tv[which(dup)[1]])))
    # values, coerced once and reported if anything was lost
    vraw <- x[[vc]]
    suppressWarnings(vnum <- as.numeric(as.character(vraw)))
    nlost <- sum(is.na(vnum) & !is.na(vraw) & nzchar(trimws(as.character(vraw))))
    if (nlost > 0)
      warning(sprintf("%d non-numeric entries in '%s' coerced to NA (decimal comma instead of point?).", nlost, vc))
    # missing combinations stay NA, so the time grid is complete and gaps are visible
    m_num <- matrix(NA_real_, length(ut), length(us), dimnames = list(NULL, us))
    m_num[cbind(ri, ci)] <- vnum
    if (isTRUE(verbose))
      message(sprintf("qc_as_xts: pivoted %d rows to %d time steps x %d stations (%.1f%% filled).",
                      nrow(x), length(ut), length(us), 100 * sum(!is.na(m_num)) / length(m_num)))
    return(xts::xts(m_num, order.by = ut))
  }
  
  #-------------------------------------------------------------------------------
  # wide layout: build the numeric matrix and hand back the xts
  
  # everything except the time column(s) is data
  m <- as.matrix(x[, setdiff(nm, tc), drop = FALSE])
  if (ncol(m) == 0) stop("No data columns left after removing the time column.")
  # report which columns are not numeric instead of coercing them into NA silently
  suppressWarnings(m_num <- matrix(as.numeric(m), nrow(m), ncol(m), dimnames = dimnames(m)))
  lost <- colnames(m)[colSums(is.na(m_num) & !is.na(m) & nzchar(trimws(m))) > 0]
  if (length(lost))
    warning(sprintf("Non-numeric entries coerced to NA in: %s (decimal comma instead of point?)",
                    paste(lost, collapse = ", ")))
  xts::xts(m_num, order.by = tv)
}
