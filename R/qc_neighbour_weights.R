#' Neighbour weight matrix of a station network (internal)
#'
#' Equations (2)-(3) of Amini et al. (2026) with the review fixes: Gaussian
#' distance weights with the bandwidth taken from the DATA-BEARING network,
#' landuse compatibility, and the k-nearest cap applied AFTER the landuse
#' masking - the published order let zero-weight stations occupy the k slots
#' and displace compatible ones. Both forest spellings count as forest.
#'
#' Validates and aligns the metadata on the way: the four required columns,
#' numeric complete coordinates, unique IDs. Stations without a metadata row
#' or without a landuse are reported with a warning, never dropped silently.
#'
#' Elevation: with a column `ELEV` (metres) in the metadata, two stations are
#' neighbours only when their elevation differs by at most `max_elev_diff`.
#' A site 400 m above the town is 2-3 K colder by lapse rate alone and would
#' be objected to against the valley floor otherwise (Goldiwil, Thun 2025).
#' A station without an elevation is not restricted; without the column the
#' rule is off and `elevation_used` says so.
#'
#' Used by `T_QC_6_spatial_consistency()` and
#' `RH_QC_7_dewpoint_consistency()`, so the network geometry exists once.
#'
#' @param metadata Data frame with ID, LAT, LON, Landuse.
#' @param stations Station names to build the network from (those present
#'   in the data); only the ones with a metadata row are used.
#' @param radius_m Maximum neighbour distance in metres.
#' @param k_neighbours At most this many nearest COMPATIBLE neighbours.
#' @param landuse_mode "graded" (same class 1, vegetated/forest pairs 0.4,
#'   else 0), "strict" (same class only) or "none" (distance only; for small
#'   networks where the landuse mask leaves nobody with neighbours).
#' @param max_elev_diff Largest elevation difference in metres between two
#'   neighbours; needs `ELEV` in the metadata.
#'
#' @return list(ids = stations used, in W's order; W = weight matrix;
#'   sigma_d = Gaussian bandwidth in metres; elevation_used = whether the
#'   elevation rule was applied). Empty ids, an empty W and sigma_d = NA when
#'   no station has a metadata row.
#'
#' @keywords internal
#' @noRd
qc_neighbour_weights <- function(metadata, stations, radius_m, k_neighbours,
                                 landuse_mode = c("graded", "strict", "none"),
                                 max_elev_diff = 100) {
  landuse_mode <- match.arg(landuse_mode)
  if (!is.numeric(max_elev_diff) || length(max_elev_diff) != 1 || max_elev_diff <= 0)
    stop("max_elev_diff must be a positive number of metres.")
  # tolerate tibbles and friends
  md <- as.data.frame(metadata)
  # the four columns the geometry relies on
  need <- c("ID", "LAT", "LON", "Landuse")
  if (!all(need %in% names(md)))
    stop(sprintf("Metadata must contain the columns %s.", paste(need, collapse = ", ")))
  # duplicated IDs would make the match ambiguous
  if (anyDuplicated(md$ID)) stop("Metadata contains duplicated IDs.")
  # only stations present in BOTH the data and the metadata can be used
  ids <- intersect(stations, md$ID)
  # stations without a metadata row can never be evaluated - say so once
  missing_md <- setdiff(stations, md$ID)
  if (length(missing_md) > 0)
    warning(sprintf("No metadata for: %s - these stations cannot be evaluated by this level.",
                    paste(missing_md, collapse = ", ")))
  # nothing matches at all: an empty network, the caller decides what that means
  if (!length(ids))
    return(list(ids = character(0), W = matrix(0, 0, 0), sigma_d = NA_real_, elevation_used = FALSE))
  # align the metadata rows to the order of the used stations
  md <- md[match(ids, md$ID), ]
  # a missing landuse would silently produce zero weights; make it visible once
  # (not in "none" mode, where the landuse is not used)
  if (landuse_mode != "none" && (anyNA(md$Landuse) || any(!nzchar(md$Landuse)))) {
    bad <- md$ID[is.na(md$Landuse) | !nzchar(md$Landuse)]
    warning(sprintf("No landuse for: %s - these stations get zero weight everywhere.",
                    paste(bad, collapse = ", ")))
    md$Landuse[is.na(md$Landuse)] <- ""
  }
  ns <- length(ids)

  # pairwise distance matrix in metres, NA diagonal (validates the coordinates)
  D <- qc_distance_matrix(md)
  # Gaussian bandwidth: median pairwise distance of the DATA-BEARING network
  sigma_d <- stats::median(D, na.rm = TRUE)
  # equation (2), distance part; the diagonal becomes weight 0
  W <- exp(-(D^2) / (2 * sigma_d^2)); W[is.na(W)] <- 0

  # the green classes - BOTH forest spellings
  veg <- c("Vegetated Areas", "Forest", "Forests")
  lu <- as.character(md$Landuse)
  # pairwise landuse factor lambda; "none" treats every pair as compatible
  LU <- if (landuse_mode == "none") matrix(1, ns, ns) else outer(lu, lu, Vectorize(function(a, b) {
    # identical non-empty class: full weight
    if (a == b && nzchar(a)) return(1)
    # graded mode only: vegetated/forest pairs get reduced weight
    if (landuse_mode == "graded" && a %in% veg && b %in% veg) return(0.4)
    # everything else: incompatible
    0
  }))
  W <- W * LU

  # the elevation rule: a pair further apart than max_elev_diff in height is
  # no pair; a station without an elevation stays unrestricted
  elevation_used <- "ELEV" %in% names(md) && any(!is.na(md$ELEV))
  if (elevation_used) {
    if (!is.numeric(md$ELEV)) stop("Metadata ELEV must be numeric (metres).")
    dz <- abs(outer(md$ELEV, md$ELEV, "-"))
    W[!is.na(dz) & dz > max_elev_diff] <- 0
  }

  # the k-nearest cap, in the FIXED order: among the compatible stations only
  for (i in seq_len(ns)) {
    # candidates = landuse-compatible AND within the radius
    cand <- which(W[i, ] > 0 & !is.na(D[i, ]) & D[i, ] <= radius_m)
    # nothing compatible in range: this station gets no support at all
    if (!length(cand)) { W[i, ] <- 0; next }
    # keep the k nearest AMONG the compatible candidates
    keep <- cand[order(D[i, cand])][seq_len(min(k_neighbours, length(cand)))]
    W[i, setdiff(seq_len(ns), keep)] <- 0
  }
  # belt and braces: nothing beyond the radius ever keeps weight
  W[!is.na(D) & D > radius_m] <- 0

  list(ids = ids, W = W, sigma_d = sigma_d, elevation_used = elevation_used)
}
