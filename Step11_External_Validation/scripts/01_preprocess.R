#!/usr/bin/env Rscript
# =============================================================================
# 01_preprocess.R
# -----------------------------------------------------------------------------
# Purpose : Preprocess E-GEOD-30999 / GSE30999 for external validation.
#   1. Read the 170 per-sample processed tables from the ArrayExpress zip.
#   2. Map probes -> gene symbols with the SAME annotation table used by the
#      discovery cohorts (Step1_DEG/probe2gene_AnnoProbe.tsv).
#   3. Collapse multiple probes per gene with the MaxMean rule: the probe with
#      the highest mean expression across all samples is retained. This is
#      identical to collapse_probes() in
#      Step1_DEG/deg_standard_pipeline_v2_AnnoProbe.R and is the rule stated
#      in the manuscript Methods (2.2.11).
#   4. Pair lesional / non-lesional biopsies by individual using the SDRF.
#
# Input  : E-GEOD-30999.processed.1.zip, E-GEOD-30999.sdrf.txt
#          (see 00_fetch_data.R; directory via GSE30999_DATA_DIR)
#          ../Step1_DEG/probe2gene_AnnoProbe.tsv
#          scripts/input/signature_18_genes.csv
# Output : data/gse30999_processed.rds        gene x sample expression matrix
#          data/gse30999_pairs.rds            81 paired samples
#          scripts/input/representative_probes.csv  18 genes -> MaxMean probes
#
# Usage  : Rscript scripts/01_preprocess.R
# =============================================================================

data_dir <- Sys.getenv("GSE30999_DATA_DIR", unset = "data")
zipfile  <- file.path(data_dir, "E-GEOD-30999.processed.1.zip")
sdrf_f   <- file.path(data_dir, "E-GEOD-30999.sdrf.txt")
stopifnot(file.exists(zipfile), file.exists(sdrf_f))
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# --- 1) SDRF: groups and complete pairs --------------------------------------
sdrf <- read.delim(sdrf_f, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
ind  <- sdrf[["Characteristics[individual]"]]
grp  <- sdrf[["FactorValue [clinical information]"]]
file <- sdrf[["Derived Array Data File"]]
stopifnot(length(ind) == length(grp), length(grp) == length(file))
cat(sprintf("[01] SDRF: %d rows | %d lesional | %d non-lesional | %d individuals\n",
            nrow(sdrf), sum(grp == "lesional"), sum(grp == "non-lesional"),
            length(unique(ind))))

by_ind   <- split(seq_along(ind), ind)
pairs_df <- do.call(rbind, lapply(by_ind, function(ix) {
  l <- file[ix][grp[ix] == "lesional"]
  n <- file[ix][grp[ix] == "non-lesional"]
  if (length(l) == 1L && length(n) == 1L)
    data.frame(patient_id = ind[ix][1],
               lesional_sample = l, nonlesional_sample = n,
               stringsAsFactors = FALSE)
}))
rownames(pairs_df) <- NULL
cat(sprintf("[01] Complete pairs: %d (manuscript target = 81)\n", nrow(pairs_df)))
if (nrow(pairs_df) != 81L)
  warning("[01] Pair count differs from 81; check the SDRF columns above.")

# --- 2) Expression matrix from the zip ---------------------------------------
cache <- file.path(tempdir(), "E-GEOD-30999_sample_tables")
if (!dir.exists(cache) || length(list.files(cache, pattern = "sample_table")) != 170L) {
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  unzip(zipfile, exdir = cache)
}
samples <- sort(unique(file))
stopifnot(length(samples) == 170L, all(file.exists(file.path(cache, samples))))

first <- read.delim(file.path(cache, samples[1]), stringsAsFactors = FALSE)
probes <- as.character(first$ID_REF)
mat <- matrix(NA_real_, nrow = length(probes), ncol = length(samples),
              dimnames = list(probes, samples))
for (f in samples) {
  d <- read.delim(file.path(cache, f), stringsAsFactors = FALSE)
  stopifnot(identical(as.character(d$ID_REF), probes))
  mat[, f] <- suppressWarnings(as.numeric(d$VALUE))
}
stopifnot(!anyNA(mat))
cat(sprintf("[01] Expression matrix: %d probes x %d samples (no NA)\n",
            nrow(mat), ncol(mat)))

# --- 3) MaxMean probe collapse (Step 1 rule) ---------------------------------
p2g_file <- file.path("..", "Step1_DEG", "probe2gene_AnnoProbe.tsv")
stopifnot(file.exists(p2g_file))
p2g <- read.delim(p2g_file, stringsAsFactors = FALSE)
p2g <- p2g[!is.na(p2g$gene) & p2g$gene != "" & p2g$probe %in% rownames(mat), , drop = FALSE]
grp_probes <- split(p2g$probe, p2g$gene)

best_probe <- function(ps) {
  if (length(ps) == 1L) return(ps)
  ps[which.max(rowMeans(mat[ps, , drop = FALSE]))]
}
sel_probe <- vapply(grp_probes, best_probe, character(1))
proc <- mat[sel_probe, , drop = FALSE]
rownames(proc) <- names(grp_probes)
proc <- proc[order(rownames(proc)), , drop = FALSE]
cat(sprintf("[01] MaxMean collapse: %d probes -> %d unique genes\n",
            nrow(p2g), nrow(proc)))

# --- 4) Record the 18 representative probes used in the reported table -------
sig   <- read.csv("scripts/input/signature_18_genes.csv", stringsAsFactors = FALSE)
genes <- sig$Gene
missing <- setdiff(genes, rownames(proc))
if (length(missing)) stop(sprintf("[01] Signature genes absent: %s",
                                  paste(missing, collapse = ", ")))
sel <- data.frame(
  Gene       = genes,
  Probe      = vapply(grp_probes[genes], best_probe, character(1)),
  n_probes   = vapply(grp_probes[genes], length, integer(1)),
  mean_expr  = round(vapply(genes, function(g)
    max(rowMeans(mat[grp_probes[[g]], , drop = FALSE])), numeric(1)), 3),
  rule       = rep("MaxMean (highest mean expression across all 170 samples)", length(genes)),
  stringsAsFactors = FALSE
)
write.csv(sel, "scripts/input/representative_probes.csv", row.names = FALSE)
cat(sprintf("[01] Wrote representative_probes.csv for %d signature genes.\n", nrow(sel)))

# --- 5) Save ------------------------------------------------------------------
saveRDS(proc,     "data/gse30999_processed.rds")
saveRDS(pairs_df, "data/gse30999_pairs.rds")
cat("[01] Saved processed matrix + pairing table to data/.\n")
cat("[01] Done. Next: Rscript scripts/02_external_validation.R\n")
