#' Temperature QC Level 6: spatial consistency (weighted neighbour consensus)
#'
#' Implements equations (1)-(4) of Amini et al. (2026): a value is removed when
#' it deviates from the Gaussian-distance- and landuse-weighted consensus of
#' its neighbours by more than max(k * sigma', delta). Three review findings
#' are corrected here:
#'
#' 1. The k-nearest cap is applied AMONG LANDUSE-COMPATIBLE stations, after
#'    the landuse masking. The published order (cap first, mask second) let
#'    zero-weight stations occupy the k slots and displace compatible ones -
#'    Lausanne Log_240 kept 1 of its 3 compatible neighbours that way, was
#'    never evaluated, and a two-week 18 K indoor episode passed unflagged.
#' 2. There is NO lower gate on sigma. The published code additionally
#'    required k*sigma > 0.5 and thereby switched the test off exactly when
#'    the neighbours agreed best (the same episode lost another 27 percent to
#'    it). The absolute floor delta alone protects against over-flagging, as
#'    in the paper's equation (4).
#' 3. Both landuse spellings "Forest" and "Forests" count as forest, because
#'    metadata and published code disagree and the mismatch silently isolates
#'    forest stations.
#'
#' sigma_d (the Gaussian bandwidth) is the median pairwise distance over the
#' stations PRESENT IN THE DATA - a deliberate, documented divergence from the
#' original, which used all metadata rows including stations without data.
#'
#' Stations that never reach the minimum number of compatible valid neighbours
#' are listed in `qc_info$t6_spatial_consistency$never_evaluated`. A station
#' this test cannot see is a coverage statement, not a clean bill.
#'
#' @param input xts of temperature, or list from a previous QC level.
#' @param metadata Data frame with columns ID, LAT, LON, Landuse. Stations are
#'   matched to data columns by ID.
#' @param k_sigma Multiplier k on the weighted neighbour spread (paper: 6).
#' @param abs_floor Absolute floor delta in Kelvin (paper: 3).
#' @param radius_m Maximum neighbour distance in metres (paper: 3000).
#' @param k_neighbours At most this many nearest COMPATIBLE neighbours (paper: 5).
#' @param min_neighbours Minimum valid neighbours to judge at all (paper: 2).
#' @param landuse_mode "graded" (paper: same class 1.0, vegetated/forest pairs
#'   0.4, else 0) or "strict" (same class only).
#' @param flag_code QC-code written by this level. Here, default is 6.
#' @param verbose Report the tally and the never-evaluated stations.
#'
#' @return The chain list with qc_data, qc_data_flagged and qc_info.
#'
#' @examples
#' \dontrun{
#' res <- T_QC_6_spatial_consistency(res, metadata = meta_biel)
#' }
#'
#' @import xts
#' @import zoo
#' @import geosphere
#' @export
T_QC_6_spatial_consistency <- function(input,
                                       metadata,
                                       k_sigma = 6,
                                       abs_floor = 3,
                                       radius_m = 3000,
                                       k_neighbours = 5,
                                       min_neighbours = 2,
                                       landuse_mode = c("graded", "strict"),
                                       flag_code = 6,
                                       verbose = TRUE) {
  #-------------------------------------------------------------------------------
  # normalise the input first and perform basic sanity checks
  
  input <- qc_prepare_input(input, what = "temperature")
  x   <- input$qc_data
  flg <- input$qc_data_flagged
  # resolve the mode choice against the two allowed values
  landuse_mode <- match.arg(landuse_mode)
  
  #-------------------------------------------------------------------------------
  # validate the parameters and the metadata, because misaligned metadata
  # would put the objections on the wrong stations
  
  # the multiplier and the floor must be positive
  if (!is.numeric(k_sigma) || k_sigma <= 0) stop("k_sigma must be a positive number.")
  if (!is.numeric(abs_floor) || abs_floor <= 0) stop("abs_floor must be a positive number of Kelvin.")
  if (!is.numeric(radius_m) || radius_m <= 0) stop("radius_m must be a positive distance in metres.")
  if (!is.numeric(min_neighbours) || min_neighbours < 1) stop("min_neighbours must be at least 1.")
  if (!is.numeric(k_neighbours) || k_neighbours < min_neighbours)
    stop("k_neighbours must be at least min_neighbours.")
  # the metadata are required for this level
  if (missing(metadata) || is.null(metadata))
    stop("This level needs the metadata (ID, LAT, LON, Landuse).")
  # tolerate tibbles and friends
  md <- as.data.frame(metadata)
  # the four columns this level relies on
  need <- c("ID", "LAT", "LON", "Landuse")
  if (!all(need %in% names(md)))
    stop(sprintf("Metadata must contain the columns %s.", paste(need, collapse = ", ")))
  # coordinates must be numeric and present
  if (!is.numeric(md$LAT) || !is.numeric(md$LON) || anyNA(md$LAT) || anyNA(md$LON))
    stop("Metadata LAT/LON must be numeric and complete.")
  # duplicated IDs would make the match ambiguous
  if (anyDuplicated(md$ID)) stop("Metadata contains duplicated IDs.")
  # only stations present in BOTH the data and the metadata can be used
  ids <- intersect(colnames(x), md$ID)
  # stations in the data without metadata can never be evaluated - say so once
  missing_md <- setdiff(colnames(x), md$ID)
  if (length(missing_md) > 0)
    warning(sprintf("No metadata for: %s - these stations cannot be evaluated by this level.",
                    paste(missing_md, collapse = ", ")))
  # without at least min_neighbours + 1 stations the condition can never be met
  if (length(ids) < min_neighbours + 1)
    stop(sprintf("Only %d station(s) match data and metadata; need at least %d.",
                 length(ids), min_neighbours + 1))
  # align the metadata rows to the order of the used stations
  md <- md[match(ids, md$ID), ]
  # a missing landuse would silently produce zero weights; make it visible once
  if (anyNA(md$Landuse) || any(!nzchar(md$Landuse))) {
    bad <- md$ID[is.na(md$Landuse) | !nzchar(md$Landuse)]
    warning(sprintf("No landuse for: %s - these stations get zero weight everywhere.",
                    paste(bad, collapse = ", ")))
    md$Landuse[is.na(md$Landuse)] <- ""
  }
  # number of stations actually used
  ns <- length(ids)
  
  #-------------------------------------------------------------------------------
  # build the static weight matrix: distance, landuse, THEN the k-nearest cap
  
  # pairwise distance matrix in metres
  D <- matrix(0, ns, ns, dimnames = list(ids, ids))
  for (i in seq_len(ns)) for (j in seq_len(ns))
    D[i, j] <- geosphere::distHaversine(c(md$LON[i], md$LAT[i]), c(md$LON[j], md$LAT[j]))
  # a station is not its own neighbour
  diag(D) <- NA
  # Gaussian bandwidth: median pairwise distance of the DATA-BEARING network
  sigma_d <- stats::median(D, na.rm = TRUE)
  # equation (2), distance part; the diagonal becomes weight 0
  W <- exp(-(D^2) / (2 * sigma_d^2)); W[is.na(W)] <- 0
  
  # the green classes - BOTH forest spellings, see header
  veg <- c("Vegetated Areas", "Forest", "Forests")
  # landuse per station as plain character
  lu <- as.character(md$Landuse)
  # pairwise landuse factor lambda
  LU <- outer(lu, lu, Vectorize(function(a, b) {
    # identical non-empty class: full weight
    if (a == b && nzchar(a)) return(1)
    # graded mode only: vegetated/forest pairs get reduced weight
    if (landuse_mode == "graded" && a %in% veg && b %in% veg) return(0.4)
    # everything else: incompatible
    0
  }))
  # combine distance and landuse
  W <- W * LU
  
  # the k-nearest cap, in the FIXED order: among the compatible stations only
  for (i in seq_len(ns)) {
    # candidates = landuse-compatible AND within the radius
    cand <- which(W[i, ] > 0 & !is.na(D[i, ]) & D[i, ] <= radius_m)
    # nothing compatible in range: this station gets no support at all
    if (!length(cand)) { W[i, ] <- 0; next }
    # keep the k nearest AMONG the compatible candidates
    keep <- cand[order(D[i, cand])][seq_len(min(k_neighbours, length(cand)))]
    # zero every slot that was not kept
    W[i, setdiff(seq_len(ns), keep)] <- 0
  }
  # belt and braces: nothing beyond the radius ever keeps weight
  W[!is.na(D) & D > radius_m] <- 0
  
  #-------------------------------------------------------------------------------
  # Perform QC Level 6
  
  # plain numeric matrix of the values (time in rows, stations in columns)
  X <- coredata(x)
  # plain numeric matrix of the flags, same shape
  previous_flag <- coredata(flg)
  # count across all time steps
  n_total <- 0
  # per-station tally for the report
  n_station <- stats::setNames(integer(ncol(X)), colnames(X))
  # track the stations this test never reached - a blind spot is not a clean bill
  never_evaluated <- stats::setNames(rep(TRUE, ns), ids)
  
  # iterate over all time steps: the weights are static, the data are not
  for (t in seq_len(nrow(X))) {
    # this time step across the used stations
    v <- X[t, ids]
    # validity mask as 0/1 for the matrix products
    vm <- as.integer(!is.na(v))
    # fewer than two values network-wide: nothing to compare
    if (sum(vm) < 2) next
    # NA-free copy so the matrix products stay defined
    vf <- ifelse(is.na(v), 0, v)
    # per station: total weight of its VALID neighbours (denominator of eq. 1)
    wsum <- as.vector(W %*% vm)
    # equation (1): weighted neighbour consensus
    mu <- as.vector(W %*% (vm * vf)) / wsum
    # a station without any weighted valid neighbour has no consensus
    mu[wsum == 0] <- NA_real_
    # per station: NUMBER of valid weighted neighbours
    nn <- rowSums((W > 0) * matrix(vm, ns, ns, byrow = TRUE))
    # who can be judged at this time step
    judged <- !is.na(v) & !is.na(mu) & nn >= min_neighbours
    if (!any(judged)) next
    # these stations have been reached by the test at least once
    never_evaluated[judged] <- FALSE
    # neighbour values as a matrix: row i holds station i's view of the network
    Xm <- matrix(rep(vf, each = ns), ns)
    # each station's consensus, repeated across its row
    Mm <- matrix(mu, ns, ns)
    # validity mask matching Xm
    vmm <- matrix(vm, ns, ns, byrow = TRUE)
    # equation (3): weighted spread of the neighbours around the consensus
    sig <- sqrt(rowSums((W * vmm) * (Xm - Mm)^2, na.rm = TRUE) / pmax(wsum, 1e-12))
    # equation (4), pure: |v - mu| > max(k*sigma, delta) - NO extra gate
    hit <- judged & (abs(v - mu) > pmax(k_sigma * sig, abs_floor))
    
    # combine the verdict with THIS time step's flag row, restricted to the
    # used stations: hit is a vector over ids, so it must meet a vector
    mask <- hit & (is.na(previous_flag[t, ids]) | previous_flag[t, ids] == 0)
    # which stations are objected to at this time step
    n_found <- sum(mask, na.rm = TRUE)
    
    # apply only if something was found
    if (n_found > 0) {
      # the station names behind the mask
      hit_ids <- ids[which(mask)]
      # blank the values so later levels never see them
      X[t, hit_ids] <- NA
      # record this level's code
      previous_flag[t, hit_ids] <- flag_code
      # add the number of new flags to the counters
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
    message(sprintf("T6 spatial consistency (max(%g*sigma, %g K), k=%d after landuse): %d flagged",
                    k_sigma, abs_floor, k_neighbours, n_total))
    if (any(never_evaluated))
      message("  never evaluated (insufficient compatible neighbours): ",
              paste(names(never_evaluated)[never_evaluated], collapse = ", "))
  }
  
  # write the updated matrices back and append this level under its own name
  input$qc_data                        <- x
  input$qc_data_flagged                <- flg
  input$qc_info$t6_spatial_consistency <- list(n_flagged            = n_total,
                                               n_flagged_by_station = n_station,
                                               sigma_d              = sigma_d,
                                               k_sigma              = k_sigma,
                                               abs_floor            = abs_floor,
                                               radius_m             = radius_m,
                                               k_neighbours         = k_neighbours,
                                               min_neighbours       = min_neighbours,
                                               landuse_mode         = landuse_mode,
                                               never_evaluated      = names(never_evaluated)[never_evaluated])
  input
}
