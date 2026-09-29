#' Humidity QC Level 1: inherit the temperature QC mask
#'
#' Temperature and relative humidity usually come from the same chip, so any
#' fault that corrupts the temperature corrupts the humidity as well. This
#' level carries the objections of the temperature run over to the humidity
#' series - cell by cell, wherever a humidity cell has a temperature flag with
#' EXACTLY the same time stamp and the same logger name. Cells without such a
#' match are left alone. The inheritance is one-directional by design.
#'
#' Without a usable temperature mask (none given, no time index, no logger
#' names, not a flag matrix, or no overlap at all) the level says so, records
#' it, and hands the pair back unchanged so the remaining levels can run.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature_flags The temperature run: its result list (preferred;
#'   its `qc_data_flagged` is taken and its `qc_info` records are merged in,
#'   so the finished humidity object documents both runs) or its flag matrix
#'   as an `xts` with logger names. NULL (default) means no inheritance.
#' @param inherit_gapfill Also inherit the plain gap-fill code 50.
#' @param verbose Report the tally and the overlap.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' t_res  <- run_qc_temperature(t_xts, metadata = meta)
#' rh_res <- RH_QC_1_inherit_temperature(rh_xts, temperature_flags = t_res)
#' rh_res <- RH_QC_1_inherit_temperature(rh_xts)   # standalone: no inheritance
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_1_inherit_temperature <- function(input,
                                        temperature_flags = NULL,
                                        inherit_gapfill = FALSE,
                                        verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks

  input <- qc_prepare_input(input, what = "humidity", level = "rh1_inherit_temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged

  #-------------------------------------------------------------------------------
  # unpack the temperature run. Anything unusable is not an error but a skip:
  # say why, record it, hand the pair back unchanged, let the chain go on

  skip <- function(reason) {
    if (isTRUE(verbose))
      message("RH1 inherit temperature mask: ", reason, " - no inheritance, continuing with the next level.")
    input$qc_info$rh1_inherit_temperature <- list(n_flagged = 0L, inherited = FALSE, reason = reason)
    input
  }
  if (is.null(temperature_flags)) return(skip("no temperature flags supplied"))
  tf <- temperature_flags
  t_info <- NULL
  if (is.list(tf) && !inherits(tf, "xts")) {
    # the whole chain result: take its mask, keep its records for the merge below
    if (!"qc_data_flagged" %in% names(tf))
      return(skip("temperature_flags list has no 'qc_data_flagged'"))
    t_info <- tf$qc_info
    tf <- tf$qc_data_flagged
  }
  # matching needs time stamps AND logger names - a bare matrix has neither
  if (!inherits(tf, "xts") || is.null(colnames(tf)))
    return(skip("temperature_flags must be an xts flag matrix with logger names"))
  # a flag matrix holds only the codes of the convention (0-9 for the levels,
  # 50-59 for the interpolation). A temperature VALUE series fails this
  M_t <- coredata(tf)
  if (!is.numeric(M_t) || !all(M_t %in% c(0:9, 50:59) | is.na(M_t)))
    return(skip("temperature_flags holds values outside the flag codes 0-9/50-59 (pass qc_data_flagged, not the values)"))

  #-------------------------------------------------------------------------------
  # the temperature codes on the humidity grid, by exact time stamp and logger
  # name (qc_match_grid): NA wherever there is no match, so those cells can
  # never inherit anything

  tm <- qc_match_grid(x, tf)
  if (is.null(tm)) return(skip("no logger name or time stamp occurs in both series"))
  loggers <- tm$loggers
  n_time_matched <- tm$n_time

  #-------------------------------------------------------------------------------
  # Perform RH QC Level 1

  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  M <- tm$M

  # every matched temperature cell that failed a test
  hit <- !is.na(M) & M > 0
  # the plain gap fill is not an objection: a missing temperature reading does
  # not make the humidity reading wrong
  if (!isTRUE(inherit_gapfill)) hit <- hit & M != 50
  # only cells that hold a humidity value and carry no earlier flag
  mask <- hit & !is.na(X) & (is.na(previous_flag) | previous_flag == 0)
  # how many cells this level objects to
  n_total <- sum(mask)
  # per-station tally for the report
  n_station <- stats::setNames(colSums(mask), colnames(X))
  # per-station coverage: humidity values that had a temperature flag to look at
  n_judged <- stats::setNames(colSums(!is.na(M) & !is.na(X)), colnames(X))
  # which temperature codes were inherited how often: the humidity flag is
  # always 1, so this is the only place the source level stays visible
  tab <- table(M[mask])
  n_by_t_code <- stats::setNames(as.integer(tab), names(tab))

  # apply only if something was found
  if (n_total > 0) {
    # blank the inherited cells so later levels never see them
    X[mask] <- NA
    # record this level's code (1 = level 1, fixed by convention)
    previous_flag[mask] <- 1
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }

  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level

  # loggers the temperature run does not know: they could not inherit anything
  unmatched <- setdiff(colnames(x), loggers)
  # report the overlap with the tally, so a low number is explainable
  if (isTRUE(verbose)) {
    message(sprintf("RH1 inherit temperature mask (gap fills %s): %d flagged, %d of %d loggers and %d of %d time steps matched",
                    if (isTRUE(inherit_gapfill)) "included" else "excluded",
                    n_total, length(loggers), ncol(X), n_time_matched, nrow(X)))
    if (length(unmatched))
      message("  no temperature for: ", paste(unmatched, collapse = ", "))
  }
  # carry the temperature run's records over, so the finished humidity object
  # documents BOTH chains; only keys that do not exist yet, nothing is overwritten
  if (!is.null(t_info))
    for (k in setdiff(names(t_info), names(input$qc_info)))
      input$qc_info[[k]] <- t_info[[k]]

  # write the updated matrices back and append this level under its own name
  input$qc_data                         <- x
  input$qc_data_flagged                 <- flg
  input$qc_info$rh1_inherit_temperature <- list(n_flagged            = n_total,
                                                n_flagged_by_station = n_station,
                                                n_judged_by_station  = n_judged,
                                                n_flagged_by_t_code  = n_by_t_code,
                                                inherited            = TRUE,
                                                inherit_gapfill      = inherit_gapfill,
                                                loggers_matched      = loggers,
                                                loggers_unmatched    = unmatched,
                                                n_time_matched       = n_time_matched)
  input
}
