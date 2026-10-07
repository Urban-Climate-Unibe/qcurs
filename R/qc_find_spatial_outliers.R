#' Weighted neighbour consensus test, time step by time step (internal)
#'
#' Equations (1), (3) and (4) of Amini et al. (2026): at every time step, each
#' station's value is compared with the weighted consensus of its valid
#' neighbours and objected to when the deviation exceeds
#' max(k_sigma * sigma', abs_floor). There is NO lower gate on sigma - the
#' published code switched the test off exactly when the neighbours agreed
#' best. A station is judged only with at least `min_neighbours` valid
#' weighted neighbours.
#'
#' Used by `T_QC_6_spatial_consistency()` (on temperature) and
#' `RH_QC_7_dewpoint_consistency()` (on the dewpoint), so the test exists once.
#'
#' @param M Numeric matrix, time in rows, stations in columns, in the order
#'   of `W`.
#' @param W Weight matrix from `qc_neighbour_weights()`.
#' @param k_sigma Multiplier on the weighted neighbour spread.
#' @param abs_floor Absolute floor of the deviation, in the unit of `M`.
#' @param min_neighbours Minimum valid weighted neighbours to judge at all.
#'
#' @return list(hit = logical matrix like `M`, TRUE at every objection;
#'   evaluated = logical per station, FALSE where the test never reached it).
#'
#' @keywords internal
#' @noRd
qc_find_spatial_outliers <- function(M, W, k_sigma, abs_floor, min_neighbours) {
  ns <- ncol(M)
  hit <- matrix(FALSE, nrow(M), ns, dimnames = dimnames(M))
  evaluated <- stats::setNames(rep(FALSE, ns), colnames(M))
  # iterate over all time steps: the weights are static, the data are not
  for (t in seq_len(nrow(M))) {
    v <- M[t, ]
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
    mu[wsum == 0] <- NA_real_
    # per station: NUMBER of valid weighted neighbours
    nn <- rowSums((W > 0) * matrix(vm, ns, ns, byrow = TRUE))
    # who can be judged at this time step
    judged <- !is.na(v) & !is.na(mu) & nn >= min_neighbours
    if (!any(judged)) next
    evaluated[judged] <- TRUE
    # equation (3): weighted spread of the neighbours around the consensus
    Xm  <- matrix(rep(vf, each = ns), ns)
    Mm  <- matrix(mu, ns, ns)
    vmm <- matrix(vm, ns, ns, byrow = TRUE)
    sig <- sqrt(rowSums((W * vmm) * (Xm - Mm)^2, na.rm = TRUE) / pmax(wsum, 1e-12))
    # equation (4), pure: |v - mu| > max(k*sigma, floor) - NO extra gate
    hit[t, ] <- judged & (abs(v - mu) > pmax(k_sigma * sig, abs_floor))
  }
  list(hit = hit, evaluated = evaluated)
}
