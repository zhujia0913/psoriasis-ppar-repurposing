#!/usr/bin/env Rscript
# =============================================================================
# 00_fetch_data.R
# -----------------------------------------------------------------------------
# Purpose : Obtain the ArrayExpress / BioStudies processed files for
#           E-GEOD-30999 (the ArrayExpress mirror of GEO series GSE30999).
#           These are the exact inputs used for the published external
#           validation statistics:
#
#             E-GEOD-30999.processed.1.zip
#                 170 per-sample tables (ID_REF / VALUE; 54,675 probes;
#                 RMA-normalized log2 values rounded to two decimals)
#             E-GEOD-30999.sdrf.txt
#                 sample/phenotype table with individual, lesion status and
#                 derived data file names
#
# REDISTRIBUTION NOTE: no source data are committed to this repository.
#   Files are read from `data/` by default, or from the directory given in
#   the environment variable GSE30999_DATA_DIR. Missing files are downloaded
#   from the EBI mirror. The .rds objects produced by 01_preprocess.R are
#   git-ignored.
#
# Usage : Rscript scripts/00_fetch_data.R
# =============================================================================

BASE <- "https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/GEOD/E-GEOD-30999"
FILES <- c("E-GEOD-30999.processed.1.zip", "E-GEOD-30999.sdrf.txt")
EXPECTED_BYTES <- c(39798039L, 109410L)
EXPECTED_MD5_SDRF <- "55d99d97c869f534a6f29ea6df85437e"

data_dir <- Sys.getenv("GSE30999_DATA_DIR", unset = "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(FILES)) {
  f <- FILES[i]
  dest <- file.path(data_dir, f)
  if (!file.exists(dest)) {
    url <- paste0(BASE, "/", f)
    cat(sprintf("[00] Downloading %s ...\n", url))
    ok <- tryCatch({
      suppressWarnings(download.file(url, dest, mode = "wb", quiet = TRUE))
      file.exists(dest)
    }, error = function(e) {
      cat("  download failed:", conditionMessage(e), "\n"); FALSE
    })
    if (!ok) {
      stop(sprintf(paste0("[00] Could not obtain %s. Download it manually from\n",
                          "     %s\n     and place it in %s/"),
                   f, BASE, data_dir))
    }
  }
  size <- file.info(dest)$size
  stopifnot(identical(as.numeric(size), as.numeric(EXPECTED_BYTES[i])))
  cat(sprintf("[00] OK  %-38s %10.0f bytes  %s\n", f, size, dest))
}

sdrf_md5 <- unname(tools::md5sum(file.path(data_dir, FILES[2])))
stopifnot(identical(sdrf_md5, EXPECTED_MD5_SDRF))
cat("[00] SDRF md5 matches the official ArrayExpress file.\n")
cat("[00] Done. Next: Rscript scripts/01_preprocess.R\n")
