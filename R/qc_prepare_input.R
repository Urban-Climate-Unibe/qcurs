#' Standard input preamble for every QC level
#'
#' Normalises the input into the canonical chain list and validates it. This is
#' the first statement of every QC level, temperature and humidity alike, so
#' that all levels share one entry contract.
#'
#' Flag convention for the whole chain, temperature and humidity alike:
#' 0 = checked and unobjected, NA = no observation, N = objected by level N,
#' 50 = plain gap filled by the optional interpolation, 50 + N = removed by
#' level N and refilled by the interpolation. The code of a level is fixed
#' by this convention and is not a parameter.
#'
#' @param input Either an `xts` object (stations in columns, one column per
#'   sensor, named after the station ID) at the start of the chain, or the
#'   `list(qc_data, qc_data_flagged, qc_info)` handed on by a previous level.
#' @param what Which chain this is: `"temperature"` or `"humidity"`. Not
#'   cosmetic - it names the records `dataset_<what>` and `conversion_<what>`,
#'   and the window based levels read their time step from there.
#'
#' @return The input as a validated chain list with three elements:
#'   \describe{
#'     \item{qc_data}{the (possibly already cleaned) data as `xts`}
#'     \item{qc_data_flagged}{the flag matrix as `xts`, same shape as `qc_data`}
#'     \item{qc_info}{the growing record; `dataset_<what>` (and on the first
#'       level `conversion_<what>`) is filled here}
#'   }
#'
#' @examples
#' \dontrun{
#' # at the start of every level:
#' input <- qc_prepare_input(input, what = "temperature")
#' x   <- input$qc_data
#' flg <- input$qc_data_flagged
#'
#' # window based levels take the resolution from here:
#' step_min <- input$qc_info$dataset_temperature$time_step_sec / 60
#' }
#'
#' @import xts
#' @import zoo
#' @export
qc_prepare_input <- function(input, what = "temperature") {

  # `what` names the records in qc_info, so only the two chain names are
  # allowed: a typo would silently create a new entry, and the window based
  # levels would not find their time step
  what <- match.arg(what, c("temperature", "humidity"))

  #-------------------------------------------------------------------------------
  # normalise the input first if necessary

  # abort if input is a data frame
  if (is.data.frame(input))
    stop("Input is a data frame. Convert it first into xts: x <- qc_as_xts(...)!")

  if (!is.list(input)) {
    # input is not a list -> this is the first level of the chain, build the pair
    if (!inherits(input, "xts"))
      stop(sprintf("Input must be an xts object or a chain result list (%s).", what))
    # never coerce silently: a character matrix would compare lexically and pass unnoticed
    if (!is.numeric(input))
      stop(sprintf("%s data must be numeric.", what))
    # refuse an empty series here, otherwise `input * 0` fails with a cryptic message
    if (nrow(input) == 0 || ncol(input) == 0)
      stop(sprintf("%s input is empty (no rows or no columns).", what))
    # take the report qc_as_xts() attached, and strip it from the data - otherwise
    # it would ride along in qc_data AND qc_data_flagged through the whole chain
    conversion <- xts::xtsAttributes(input)$qc_conversion
    if (!is.null(conversion)) xts::xtsAttributes(input) <- NULL
    # build the list
    input <- list(qc_data = input,
                  # generates a 0 matrix (same dimension as qc_data). All entries are 0 except NA
                  qc_data_flagged = input * 0,
                  qc_info = list())
    # file the report under the chain's name, next to dataset_<what>
    if (!is.null(conversion)) input$qc_info[[paste0("conversion_", what)]] <- conversion
  }

  #-------------------------------------------------------------------------------
  # then validate the input. From here on input is always a list: either it came
  # as one, or the block above just built it

  # ...check whether the list contains correctly labelled data
  if (!all(c("qc_data", "qc_data_flagged") %in% names(input)))
    stop("Input list must contain 'qc_data' and 'qc_data_flagged'.")
  # if yes, then assign the data to its corresponding variable
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  # ... check whether qc_data is an xts object
  if (!inherits(x, "xts"))
    stop("qc_data must be an xts object.")
  # ... check whether qc_data_flagged is an xts object
  if (!inherits(flg, "xts"))
    stop("qc_data_flagged must be an xts object.")
  # ... check whether qc_data is numeric
  if (!is.numeric(x))
    stop("qc_data must be numeric.")
  # ... check whether qc_data_flagged is numeric
  if (!is.numeric(flg))
    stop("qc_data_flagged must be numeric.")
  # ... check whether there is any data at all (empty loops would report "0 flagged" silently)
  if (nrow(x) == 0 || ncol(x) == 0)
    stop("qc_data is empty (no rows or no columns).")
  # ... check whether qc_data and qc_data_flagged have the same dimension
  if (!identical(dim(x), dim(flg)))
    stop("qc_data and qc_data_flagged must have identical dimensions.")
  # ... check whether the stations are named at all (unnamed columns make every station loop a no-op)
  if (is.null(colnames(x)) || any(is.na(colnames(x))) || any(!nzchar(colnames(x))))
    stop("qc_data needs non-empty column names (station IDs).")
  # ... check whether every station name occurs only once (duplicates make x[, name] ambiguous)
  if (anyDuplicated(colnames(x))) stop("qc_data has duplicated station names.")
  # ... check whether qc_data and qc_data_flagged have the same column-names
  if (!identical(colnames(x), colnames(flg)))
    stop("qc_data and qc_data_flagged must have identical station columns in identical order.")
  # ... check whether qc_data and qc_data_flagged have the same time-index
  if (!isTRUE(all.equal(as.numeric(index(x)), as.numeric(index(flg)))))
    stop("qc_data and qc_data_flagged must have identical time indices.")
  # ... check whether any time stamp occurs twice (breaks rolling windows and time based subsetting)
  if (anyDuplicated(index(x)))
    stop("qc_data has duplicated time stamps.")
  # ... check whether the time index runs forward (rolling windows assume ascending order)
  if (is.unsorted(index(x)))
    stop("qc_data time index must be in ascending order.")

  # ... determine the dominant time step ONCE, so window based levels do not have to
  #     guess it from the first two stamps, and warn if the series is not on a regular grid
  time_step_sec <- NA_real_
  irregular <- FALSE
  # the dataset entry is named by the variable, so a temperature run and a
  # humidity run on the SAME chain list keep separate records:
  # qc_info$dataset_temperature and qc_info$dataset_humidity
  ds_name <- paste0("dataset_", what)
  # the grid is a property of the DATASET, not of the level: warn on the first
  # level of EACH variable only, otherwise the same warning appears once per
  # level and drowns the real ones
  first_pass <- is.null(input$qc_info[[ds_name]])
  if (nrow(x) > 1) {
    # spacing between consecutive time stamps, in seconds
    steps <- diff(as.numeric(index(x)))
    # the most frequent spacing is the nominal resolution of the network
    tab <- table(steps)
    time_step_sec <- as.numeric(names(tab)[which.max(tab)])
    # more than one distinct spacing means gaps or a mixed resolution
    irregular <- length(tab) > 1
    if (irregular && first_pass)
      warning(sprintf("Irregular time steps: %d of %d intervals differ from the dominant %g min.",
                      sum(steps != time_step_sec), length(steps), time_step_sec / 60))
  }

  # ... check whether qc_info exists. if not, create it
  if (is.null(input$qc_info)) input$qc_info <- list()
  # ... record what we learned about the dataset, so every level can rely on it.
  #     ONE entry per VARIABLE, rewritten on every call of that variable's
  #     chain: within a chain the dataset never changes, so the rewrite is
  #     idempotent - while a later humidity run writes its own entry next to
  #     the temperature one instead of overwriting it.
  input$qc_info[[ds_name]] <- list(variable             = what,
                                   n_time               = nrow(x),
                                   n_stations           = ncol(x),
                                   stations             = colnames(x),
                                   time_start           = index(x)[1],
                                   time_end             = index(x)[nrow(x)],
                                   time_step_sec        = time_step_sec,
                                   irregular_time_steps = irregular)

  # hand the validated chain list back to the calling level
  input
}
