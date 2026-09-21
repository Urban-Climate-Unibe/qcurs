#' Humidity QC Level 7: dewpoint consistency with the neighbouring stations
#'
#' The spatial comparison runs on the DEWPOINT, not on relative humidity: RH
#' is spatially noisy because it carries the temperature field with it, while
#' the dewpoint is a smooth tracer of the air mass. Dewpoint after Magnus with
#' the Alduchov & Eskridge (1996) coefficients. The weighting mirrors
#' temperature level 6 including all its fixes: k-nearest cap AFTER the
#' landuse masking, NO sigma gate, both forest spellings. A "Td <= T" check is
#' deliberately absent: with Td derived from T and RH it holds automatically
#' for RH <= 100, and RH above 100 is caught by level 2.
#'
#' Objections land on the HUMIDITY cells: the dewpoint is only the yardstick,
#' the humidity is the measured quantity under test. Stations that never reach
#' the minimum number of compatible valid neighbours are listed in
#' `qc_info$rh7_dewpoint_consistency$never_evaluated`.
#'
#' @param input xts of relative humidity, or list from a previous QC level.
#' @param temperature The station temperatures: an xts with the same time
#'   index, or the whole result list of the temperature chain (its cleaned
#'   `qc_data` is taken). Stations are matched by name.
#' @param metadata Data frame with columns ID, LAT, LON, Landuse.
#' @param td_k Multiplier on the weighted neighbour spread of the dewpoint.
#' @param td_floor Absolute floor in Kelvin of dewpoint deviation.
#' @param td_radius_m Maximum neighbour distance in metres.
#' @param td_k_neighbours At most this many nearest COMPATIBLE neighbours.
#' @param td_min_neighbours Minimum valid neighbours to judge at all.
#' @param verbose Report the tally and the never-evaluated stations.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- RH_QC_7_dewpoint_consistency(res, temperature = t_res, metadata = meta)
#' }
#'
#' @import xts
#' @import zoo
#' @import geosphere
#' @export
RH_QC_7_dewpoint_consistency <- function(input,
                                         temperature,
                                         metadata,
                                         td_k = 6,
                                         td_floor = 3,
                                         td_radius_m = 3000,
                                         td_k_neighbours = 5,
                                         td_min_neighbours = 2,
                                         verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "humidity")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  
  #-------------------------------------------------------------------------------
  # validate the parameters, the temperature series and the metadata
  
  if (!is.numeric(td_k) || td_k <= 0) stop("td_k must be a positive number.")
  if (!is.numeric(td_floor) || td_floor <= 0) stop("td_floor must be a positive number of Kelvin.")
  if (!is.numeric(td_radius_m) || td_radius_m <= 0) stop("td_radius_m must be a positive distance in metres.")
  if (!is.numeric(td_min_neighbours) || td_min_neighbours < 1) stop("td_min_neighbours must be at least 1.")
  if (!is.numeric(td_k_neighbours) || td_k_neighbours < td_min_neighbours)
    stop("td_k_neighbours must be at least td_min_neighbours.")
  # the temperature series is required to compute the dewpoint
  if (missing(temperature) || is.null(temperature)) stop("This level needs the temperature series.")
  tx <- temperature
  if (is.list(tx) && !is.data.frame(tx) && !inherits(tx, "xts")) {
    if (!"qc_data" %in% names(tx)) stop("temperature list must contain 'qc_data'.")
    tx <- tx$qc_data
  }
  if (is.null(dim(tx))) stop("temperature must be an xts, a matrix, or a chain result list.")
  if (nrow(tx) != nrow(x))
    stop(sprintf("Row mismatch: humidity has %d time steps, temperature has %d.", nrow(x), nrow(tx)))
  if (inherits(tx, "xts") &&
      !isTRUE(all.equal(as.numeric(zoo::index(tx)), as.numeric(zoo::index(x)))))
    stop("Time index of the temperature series differs from the humidity series.")
  if (is.null(colnames(tx))) stop("temperature needs station names to be matched to the humidity columns.")
  # the metadata are required for the neighbour geometry
  if (missing(metadata) || is.null(metadata)) stop("This level needs the metadata (ID, LAT, LON, Landuse).")
  md <- as.data.frame(metadata)
  need <- c("ID", "LAT", "LON", "Landuse")
  if (!all(need %in% names(md)))
    stop(sprintf("Metadata must contain the columns %s.", paste(need, collapse = ", ")))
  if (anyDuplicated(md$ID)) stop("Metadata contains duplicated IDs.")
  # usable stations need humidity, temperature AND metadata
  ids <- intersect(intersect(colnames(x), colnames(tx)), md$ID)
  # a consensus over fewer than 3 stations is not a network statement:
  # refuse audibly and hand the pair back UNCHANGED, in the standard shape
  if (length(ids) < 3) {
    if (isTRUE(verbose))
      message("RH7 dewpoint consistency skipped: fewer than 3 stations matched humidity, temperature and metadata.")
    input$qc_info$rh7_dewpoint_consistency <- list(n_flagged = 0L, skipped = TRUE,
                                                   reason = "fewer than 3 matched stations")
    return(input)
  }
  md <- md[match(ids, md$ID), ]
  # a missing landuse silently produces zero weights; make it visible once
  if (anyNA(md$Landuse) || any(!nzchar(md$Landuse))) {
    bad <- md$ID[is.na(md$Landuse) | !nzchar(md$Landuse)]
    warning(sprintf("No landuse for: %s - these stations get zero weight everywhere.",
                    paste(bad, collapse = ", ")))
    md$Landuse[is.na(md$Landuse)] <- ""
  }
  ns <- length(ids)
  
  #-------------------------------------------------------------------------------
  # build the static weight matrix: distance, landuse, THEN the k-nearest cap
  
  # pairwise distance matrix in metres
  D <- matrix(0, ns, ns, dimnames = list(ids, ids))
  for (i in seq_len(ns)) for (j in seq_len(ns))
    D[i, j] <- geosphere::distHaversine(c(md$LON[i], md$LAT[i]), c(md$LON[j], md$LAT[j]))
  # a station is not its own neighbour
  diag(D) <- NA
  # Gaussian bandwidth from the used network
  sigma_d <- stats::median(D, na.rm = TRUE)
  # distance weights; the diagonal becomes 0
  W <- exp(-(D^2) / (2 * sigma_d^2)); W[is.na(W)] <- 0
  # the green classes - BOTH forest spellings
  veg <- c("Vegetated Areas", "Forest", "Forests")
  lu <- as.character(md$Landuse)
  # pairwise landuse factor
  LU <- outer(lu, lu, Vectorize(function(a, b) {
    if (a == b && nzchar(a)) return(1)
    if (a %in% veg && b %in% veg) return(0.4)
    0
  }))
  W <- W * LU
  # the k-nearest cap, in the FIXED order: among the compatible stations only
  for (i in seq_len(ns)) {
    cand <- which(W[i, ] > 0 & !is.na(D[i, ]) & D[i, ] <= td_radius_m)
    if (!length(cand)) { W[i, ] <- 0; next }
    keep <- cand[order(D[i, cand])][seq_len(min(td_k_neighbours, length(cand)))]
    W[i, setdiff(seq_len(ns), keep)] <- 0
  }
  W[!is.na(D) & D > td_radius_m] <- 0
  
  #-------------------------------------------------------------------------------
  # the dewpoint field the test runs on
  
  # Magnus formula with the Alduchov & Eskridge (1996) coefficients
  dewpoint <- function(t_c, rh_pct) {
    a <- 17.625; b <- 243.04
    # clamp into the defined domain of the logarithm
    r <- pmin(pmax(rh_pct, 0.1), 100)
    g <- log(r / 100) + (a * t_c) / (b + t_c)
    (b * g) / (a - g)
  }
  
  #-------------------------------------------------------------------------------
  # Perform RH QC Level 7
  
  # plain numeric matrix of the humidity values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # temperature matrix aligned to the used stations
  TT <- as.matrix(coredata(tx))[, ids, drop = FALSE]; storage.mode(TT) <- "double"
  # dewpoint per cell, from the CURRENT (cleaned) humidity and temperature
  TD <- dewpoint(TT, X[, ids, drop = FALSE])
  # count across all time steps
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # track the stations this test never reached
  never_evaluated <- stats::setNames(rep(TRUE, ns), ids)
  
  # iterate over all time steps: the weights are static, the data are not
  for (t in seq_len(nrow(X))) {
    # the dewpoints at this time step
    v <- TD[t, ]
    # validity mask as 0/1 for the matrix products
    vm <- as.integer(!is.na(v))
    # fewer than two values network-wide: nothing to compare
    if (sum(vm) < 2) next
    # NA-free copy so the matrix products stay defined
    vf <- ifelse(is.na(v), 0, v)
    # per station: total weight of its VALID neighbours
    wsum <- as.vector(W %*% vm)
    # the weighted neighbour consensus of the dewpoint
    mu <- as.vector(W %*% (vm * vf)) / wsum
    mu[wsum == 0] <- NA_real_
    # per station: NUMBER of valid weighted neighbours
    nn <- rowSums((W > 0) * matrix(vm, ns, ns, byrow = TRUE))
    # who can be judged at this time step
    judged <- !is.na(v) & !is.na(mu) & nn >= td_min_neighbours
    if (!any(judged)) next
    never_evaluated[judged] <- FALSE
    # the weighted spread of the neighbours around the consensus
    Xm  <- matrix(rep(vf, each = ns), ns)
    Mm  <- matrix(mu, ns, ns)
    vmm <- matrix(vm, ns, ns, byrow = TRUE)
    sig <- sqrt(rowSums((W * vmm) * (Xm - Mm)^2, na.rm = TRUE) / pmax(wsum, 1e-12))
    # objection when the dewpoint deviation exceeds max(k*sigma, floor) - NO gate
    hit <- judged & (abs(v - mu) > pmax(td_k * sig, td_floor))
    
    # combine the verdict with THIS time step's flag row over the used stations
    mask <- hit & (is.na(previous_flag[t, ids]) | previous_flag[t, ids] == 0)
    # how many stations are objected to at this time step
    n_found <- sum(mask, na.rm = TRUE)
    
    # apply only if something was found - onto the HUMIDITY, not the dewpoint
    if (n_found > 0) {
      hit_ids <- ids[which(mask)]
      X[t, hit_ids] <- NA
      previous_flag[t, hit_ids] <- 7
      n_station[hit_ids] <- n_station[hit_ids] + 1L
      n_total <- n_total + n_found
    }
  }
  
  # write the matrices back into the xts shells, keeping index and column names
  if (n_total > 0) {
    x[]   <- X
    flg[] <- previous_flag
  }
  
  #-------------------------------------------------------------------------------
  # report and hand the pair on to the next level
  
  # report so a zero-hit run is visibly a run, not a skip - and name the blind spots
  if (isTRUE(verbose)) {
    message(sprintf("RH7 dewpoint consistency (max(%g*sigma, %g K), k=%d after landuse): %d flagged",
                    td_k, td_floor, td_k_neighbours, n_total))
    if (any(never_evaluated))
      message("  never evaluated (insufficient compatible neighbours): ",
              paste(names(never_evaluated)[never_evaluated], collapse = ", "))
  }
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                          <- x
  input$qc_data_flagged                  <- flg
  input$qc_info$rh7_dewpoint_consistency <- list(n_flagged            = n_total,
                                                 n_flagged_by_station = n_station,
                                                 sigma_d              = sigma_d,
                                                 td_k                 = td_k,
                                                 td_floor             = td_floor,
                                                 td_radius_m          = td_radius_m,
                                                 td_k_neighbours      = td_k_neighbours,
                                                 td_min_neighbours    = td_min_neighbours,
                                                 never_evaluated      = names(never_evaluated)[never_evaluated])
  input
}
