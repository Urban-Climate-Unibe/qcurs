#' Summarise a QC result per level and station
#'
#' One row per level and station: how many values the level objected to, how
#' many it actually reached, and that reach as a share of the station's
#' observed values. `judged_pct` of 0 means the level never saw the station -
#' a flag of 0 there is no clean bill. A level that skipped itself leaves no
#' rows; its reason is in `qc_info`. The interpolation leaves none either, it
#' adds values rather than judging them. Records of the other chain that
#' travelled along (the temperature run inside a humidity result) are left out.
#'
#' @param res A chain result list (qc_data, qc_data_flagged, qc_info).
#' @param chain `"t"` or `"rh"`; NULL picks the chain the result belongs to.
#'
#' @return Data frame with `level`, `station`, `flagged`, `judged`,
#'   `judged_pct`, in the order the levels ran.
#'
#' @examples
#' \dontrun{
#' s <- qc_summary(t_res)
#' xtabs(flagged ~ station + level, s)
#' xtabs(judged_pct ~ station + level, s)
#' }
#'
#' @import zoo
#' @export
qc_summary <- function(res, chain = NULL) {
  # only a chain result carries the three parts this table is built from
  if (!is.list(res) || !all(c("qc_data", "qc_data_flagged", "qc_info") %in% names(res)))
    stop("res must be a chain result list with qc_data, qc_data_flagged and qc_info.")
  # the chain this result belongs to: humidity once a humidity level touched it
  if (is.null(chain)) chain <- if (!is.null(res$qc_info$dataset_humidity)) "rh" else "t"
  chain <- match.arg(chain, c("t", "rh"))
  # the level records of that chain, in the order they ran
  keys <- grep(paste0("^", chain, "[0-9]+_"), names(res$qc_info), value = TRUE)
  rows <- lapply(keys, function(k) {
    r <- res$qc_info[[k]]
    # skipped levels and the interpolation have no per-station record
    if (is.null(r$n_flagged_by_station)) return(NULL)
    data.frame(level   = k,
               station = names(r$n_flagged_by_station),
               flagged = as.integer(r$n_flagged_by_station),
               judged  = as.integer(r$n_judged_by_station),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # nothing ran: an empty table with the right columns
  if (is.null(out))
    return(data.frame(level = character(0), station = character(0), flagged = integer(0),
                      judged = integer(0), judged_pct = numeric(0), stringsAsFactors = FALSE))
  # observed values per station: the flag matrix is NA where nothing was
  # measured and 50 where a plain gap was filled; everything else was observed
  F <- zoo::coredata(res$qc_data_flagged)
  n_obs <- colSums(!is.na(F) & F != 50)
  out$judged_pct <- round(100 * out$judged / pmax(n_obs[out$station], 1), 1)
  rownames(out) <- NULL
  out
}
