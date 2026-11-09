#!/usr/bin/env Rscript
# =============================================================================
# 00_fetch_data.R
# -----------------------------------------------------------------------------
# Purpose : Obtain the EXTERNAL VALIDATION cohort GSE30999 from NCBI GEO.
#
# Dataset : GSE30999
#           - 85 patients with moderate-to-severe psoriasis
#           - 170 paired lesional / non-lesional skin biopsies
#           - Affymetrix HG-U133 Plus 2.0 (GPL570)
#           - After QC, 81 complete pairs were retained for analysis
#           (see manuscript Methods 2.2.11 / Results 3.11)
#
# REDISTRIBUTION NOTE:
#   Raw/processed data from GEO are NOT redistributed inside this repository.
#   Run this script (requires internet access + Bioconductor) to fetch them.
#   The resulting .rds files are git-ignored (see .gitignore at repo root).
#
# Output  : data/gse30999_expression.rds  - expression matrix, probes x samples (log2, RMA)
#           data/gse30999_pheno.rds       - data.frame of sample metadata
#
# Usage   : Rscript scripts/00_fetch_data.R
#
# Provenance: external-validation stage of the psoriasis PPAR drug-repurposing
#             study. Journal-independent: no target journal is assumed.
#             Verified against R 4.5.0 / GEOquery 2.80.0; getGSEMatrix() was
#             removed from current GEOquery, so getGEO(GSEMatrix = TRUE) is used.
# =============================================================================

suppressMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")
  library(GEOquery)
})

gse_id  <- "GSE30999"
destdir <- "data"
dir.create(destdir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[00] Downloading %s from GEO (this may take a few minutes)...\n", gse_id))
local_gz <- file.path(destdir, paste0(gse_id, "_series_matrix.txt.gz"))
if (file.exists(local_gz)) {
  cat(sprintf("[00] Using local series-matrix file: %s\n", local_gz))
  gse_list <- getGEO(filename = local_gz, GSEMatrix = TRUE, getGPL = FALSE)
} else {
  gse_list <- getGEO(GEO = gse_id, GSEMatrix = TRUE, getGPL = FALSE,
                     destdir = destdir)
}
gse <- if (is.list(gse_list) && !inherits(gse_list, "ExpressionSet")) gse_list[[1]] else gse_list

# --- Extract expression matrix and pheno data -------------------------------
# getGEO(GSEMatrix = TRUE) returns a list of ExpressionSet objects (one per platform).
if (inherits(gse, "ExpressionSet")) {
  expr.mat <- Biobase::exprs(gse)      # probes x samples
  pdat     <- Biobase::pData(gse)      # sample metadata
} else if (is.list(gse) && !is.data.frame(gse)) {
  mat_df   <- gse[[1]]                 # series-matrix table (first column = probe)
  expr.mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(expr.mat) <- mat_df[[1]]
  pdat     <- data.frame(sample_id = colnames(expr.mat), stringsAsFactors = FALSE)
  # enrich with any phenoData attribute if present
  if (!is.null(attr(mat_df, "phenoData")))
    pdat <- cbind(pdat, attr(mat_df, "phenoData"))
} else {
  stop("[00] Unexpected object returned by getGEO(); inspect `gse` manually.")
}

cat(sprintf("[00] Retrieved %d probes x %d samples.\n", nrow(expr.mat), ncol(expr.mat)))
cat("[00] Sample metadata columns:\n")
print(colnames(pdat))

saveRDS(expr.mat, file.path(destdir, "gse30999_expression.rds"))
saveRDS(pdat,    file.path(destdir, "gse30999_pheno.rds"))
cat(sprintf("[00] Saved expression + pheno to %s/.\n", destdir))

# -----------------------------------------------------------------------------
# OPTIONAL: re-derive RMA from raw CEL files (instead of the submitter matrix).
# Uncomment if you prefer to normalize from scratch (matches "RMA-normalized"
# wording in the manuscript). Requires the affy package and the .CEL.gz files
# placed under data/CEL/.
#
#   if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
#   library(affy)
#   ab       <- ReadAffy(celfile.path = file.path(destdir, "CEL"))
#   expr.mat <- exprs(rma(ab))   # RMA, log2 scale
#   saveRDS(expr.mat, file.path(destdir, "gse30999_expression.rds"))
# -----------------------------------------------------------------------------
