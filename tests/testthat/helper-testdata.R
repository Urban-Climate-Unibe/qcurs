# Shared fixtures: the Biel test dataset shipped in inst/extdata, loaded once
# per test run. The tidy file is pivoted here so every test file works on the
# same wide xts the chain expects.

biel_tidy <- function() {
  # ASCII-only file, so a plain read is safe here (see biel_meta for the locale note)
  utils::read.csv(system.file("extdata", "Biel_test_tidy.csv", package = "qcurbnet"),
                  sep = ";", dec = ".", stringsAsFactors = FALSE)
}

biel_meta <- function() {
  # Read via readLines + iconv rather than a UTF-8 connection: read.csv() on a
  # UTF-8 connection dies with "invalid input" under a C locale, which is what
  # R CMD check and many CI runners use. Bytes first, encoding second.
  path <- system.file("extdata", "metadata.csv", package = "qcurbnet")
  txt  <- readLines(path, warn = FALSE)
  txt  <- iconv(txt, from = "UTF-8", to = "UTF-8", sub = "")   # drop undecodable bytes
  txt[1] <- sub("^\ufeff", "", txt[1])                         # strip the BOM from the header
  md <- utils::read.csv(text = txt, sep = ";", dec = ".", stringsAsFactors = FALSE)
  names(md)[1] <- "ID"
  md[md$City == "BIEL", ]
}

biel_manifest <- function() {
  utils::read.csv(system.file("extdata", "planted_errors.csv", package = "qcurbnet"),
                  sep = ";", stringsAsFactors = FALSE)
}

# the full chain, run once and memoised, because it takes a few seconds
.chain_cache <- new.env(parent = emptyenv())
biel_chain <- function() {
  if (is.null(.chain_cache$res)) {
    x <- qc_as_xts(biel_tidy(), verbose = FALSE)
    md <- biel_meta()
    res <- T_QC_1_gross_error(x, verbose = FALSE)
    res <- T_QC_2_out_of_range(res, verbose = FALSE)
    res <- T_QC_3_time_consistency(res, verbose = FALSE)
    res <- T_QC_4_temporal_persistence(res, verbose = FALSE)
    res <- T_QC_5_climatic_outliers(res, verbose = FALSE)
    res <- T_QC_6_spatial_consistency(res, metadata = md, verbose = FALSE)
    res <- T_QC_7_spatiotemporal_consistency(res, metadata = md, verbose = FALSE)
    res <- T_QC_8_diurnal_range(res, verbose = FALSE)
    res <- T_QC_9_interpolate(res, verbose = FALSE)
    .chain_cache$x <- x
    .chain_cache$res <- res
  }
  .chain_cache
}

# manifest rows -> cell indices in the flag matrix, keyed by (time, station)
planted_cells <- function(id, flg) {
  m <- biel_manifest(); m <- m[m$id == id, ]
  ri <- match(m$time, format(zoo::index(flg), "%Y-%m-%d %H:%M:%S"))
  ci <- match(m$station, colnames(flg))
  testthat::expect_false(anyNA(ri) || anyNA(ci))
  cbind(ri, ci)
}
