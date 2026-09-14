#' Humidity QC Level 1: inherit the temperature QC mask
#'
#' Temperature and relative humidity usually come from the same chip, so any
#' fault that corrupts the temperature corrupts the humidity as well. This
#' level therefore carries the objections of the temperature run over to the
#' humidity series. The inheritance is one-directional by design: a humidity
#' flag is never pushed back onto the temperature.
#'
#' Which temperature codes are inherited: every code above 0 counts as "this
#' cell failed a test", including the interpolation REFILL codes 51-58, which
#' mark cells a QC level had removed. The plain gap-fill code 50 is NOT
#' inherited by default - a missing temperature reading is not evidence
#' against the humidity reading. Set `inherit_gapfill = TRUE` to include it.
#'
#' The mask must describe THE SAME dataset: identical number of time steps,
#' identical time index, identical station set. Anything else aborts, because
#' a silently misaligned mask flags the wrong cells and nothing looks wrong
#' afterwards. Column ORDER may differ - it is aligned by name.
#'
#' Like every level, this one runs standalone. Called WITHOUT
#' `temperature_flags` it does nothing, says so, and records the refusal in
#' `qc_info` - so a humidity chain can be run before its temperature
#' counterpart exists without the result silently pretending inheritance had
#' happened.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature_flags The temperature run's flag matrix: an `xts`, a
#'   plain matrix, or the whole result list of the temperature chain - then
#'   its `qc_data_flagged` is taken AND its `qc_info` records are merged into
#'   this chain, so the finished humidity object documents both runs. NULL
#'   (default) means no inheritance.
#' @param inherit_gapfill Also inherit the plain gap-fill code 50.
#' @param flag_code QC-code written by this level. Here, default is 1.
#' @param verbose Report the tally.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' t_res  <- T_QC_8_diurnal_range(t_res)          # the finished temperature run
#' rh_res <- RH_QC_1_inherit_temperature(rh_xts, temperature_flags = t_res)
#'
#' rh_res <- RH_QC_1_inherit_temperature(rh_xts)  # standalone: no inheritance
#' }
#'
#' @import xts
#' @import zoo
#' @export
RH_QC_1_inherit_temperature <- function(input,
                                        temperature_flags = NULL,
                                        inherit_gapfill = FALSE,
                                        flag_code = 1,
                                        verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # without a temperature mask there is nothing to inherit: say so and hand the
  # pair back UNCHANGED, so this level can also be run on its own
  
  if (is.null(temperature_flags)) {
    if (isTRUE(verbose))
      message("RH1 inherit temperature mask: no temperature flags supplied - no inheritance possible.")
    # record the refusal, so a zero-flag run is distinguishable from a clean one
    input$qc_info$rh1_inherit_temperature <- list(n_flagged       = 0L,
                                                  inherited       = FALSE,
                                                  reason          = "no temperature_flags supplied",
                                                  inherit_gapfill = inherit_gapfill)
    return(input)
  }
  
  #-------------------------------------------------------------------------------
  # unpack the mask argument: xts, plain matrix, or a whole chain result list
  
  tf <- temperature_flags
  if (is.list(tf) && !is.data.frame(tf) && !inherits(tf, "xts")) {
    # a full chain result was passed: it must actually contain the mask
    if (!"qc_data_flagged" %in% names(tf))
      stop("temperature_flags list must contain 'qc_data_flagged'.")
    # carry the temperature run's records over, so the finished humidity object
    # documents BOTH chains (dataset_temperature, t1_..., dataset_humidity,
    # rh1_...). Only keys that do not exist yet are copied - nothing of the
    # humidity chain is ever overwritten.
    if (!is.null(tf$qc_info))
      for (k in setdiff(names(tf$qc_info), names(input$qc_info)))
        input$qc_info[[k]] <- tf$qc_info[[k]]
    tf <- tf$qc_data_flagged
  }
  # anything else must at least be rectangular and numeric
  if (is.null(dim(tf))) stop("temperature_flags must be an xts, a matrix, or a chain result list.")
  
  #-------------------------------------------------------------------------------
  # strict alignment: the same dataset or nothing
  
  # the number of time steps must match exactly
  if (nrow(tf) != nrow(x))
    stop(sprintf("Row mismatch: humidity has %d time steps, temperature mask has %d.",
                 nrow(x), nrow(tf)))
  # the number of stations must match exactly
  if (ncol(tf) != ncol(x))
    stop(sprintf("Column mismatch: humidity has %d stations, temperature mask has %d.",
                 ncol(x), ncol(tf)))
  # if the mask carries a time index, it must be the same one
  if (inherits(tf, "xts") &&
      !isTRUE(all.equal(as.numeric(zoo::index(tf)), as.numeric(zoo::index(x)))))
    stop("Time index of the temperature mask differs from the humidity series.")
  # if the mask carries station names, the SETS must be identical ...
  if (!is.null(colnames(tf))) {
    if (!setequal(colnames(tf), colnames(x)))
      stop(sprintf("Station names differ. Only in humidity: %s | only in mask: %s",
                   paste(setdiff(colnames(x), colnames(tf)), collapse = ", "),
                   paste(setdiff(colnames(tf), colnames(x)), collapse = ", ")))
    # ... the ORDER may differ, so align it by name
    tf <- tf[, colnames(x)]
  } else {
    # no names means we cannot verify the assignment - warn, then trust the order
    warning("temperature_flags has no station names; inheriting by column POSITION.")
  }
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 1
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # plain numeric matrix of the TEMPERATURE codes, same shape
  M <- as.matrix(coredata(tf)); storage.mode(M) <- "double"
  
  # every temperature cell that failed a test. This is elementwise over the
  # whole matrix - unlike the temperature levels there is no per-station
  # computation here, so no station loop is needed; the mask line below has
  # the same shape, only both operands are matrices instead of vectors.
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
  
  # apply only if something was found
  if (n_total > 0) {
    # blank the inherited cells so later levels never see them
    X[mask] <- NA
    # record this level's code
    previous_flag[mask] <- flag_code
    # write the matrices back into the xts shells, keeping index and column names
    x[]   <- X
    flg[] <- previous_flag
  }
  
  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level
  
  # report so a zero-hit run is visibly a run, not a skip
  if (isTRUE(verbose))
    message(sprintf("RH1 inherit temperature mask (gap fills %s): %d flagged",
                    if (isTRUE(inherit_gapfill)) "included" else "excluded", n_total))
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                          <- x
  input$qc_data_flagged                  <- flg
  input$qc_info$rh1_inherit_temperature  <- list(n_flagged            = n_total,
                                                 n_flagged_by_station = n_station,
                                                 inherited            = TRUE,
                                                 inherit_gapfill      = inherit_gapfill)
  input
}
