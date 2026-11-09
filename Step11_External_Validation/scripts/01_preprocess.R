#!/usr/bin/env Rscript
# =============================================================================
# 01_preprocess.R
# -----------------------------------------------------------------------------
# Purpose : Preprocess GSE30999 for external validation.
#   1. Map Affymetrix HG-U133 Plus 2.0 probes -> gene symbols.
#   2. Collapse multiple probes mapping to the same gene to the MEDIAN
#      (manuscript 2.2.11: "multiple probes mapping to the same gene
#       collapsed to the median").
#   3. Parse sample metadata to identify lesional vs non-lesional status and
#      the patient-level pairing; keep only complete pairs (target: 81 pairs).
#
# Input   : data/gse30999_expression.rds, data/gse30999_pheno.rds
#           (produced by 00_fetch_data.R)
# Output  : data/gse30999_processed.rds  - gene-symbol x sample expression matrix
#           data/gse30999_pairs.rds      - pairing table (patient/lesional/non-lesional)
#
# Usage   : Rscript scripts/01_preprocess.R
#
# NOTE: the lesion-status and patient-ID parsing below is dataset-specific.
#       Run 00_fetch_data.R first and inspect the metadata columns printed by
#       it in data/gse30999_pheno.rds before trusting the pairing. The script
#       warns if the number of complete pairs differs from the reported 81.
# =============================================================================

suppressMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("hgu133plus2.db", quietly = TRUE))
    BiocManager::install("hgu133plus2.db")
  library(hgu133plus2.db)
})

expr.mat <- readRDS("data/gse30999_expression.rds")
pdat     <- readRDS("data/gse30999_pheno.rds")

# --- 1) Probe -> gene symbol annotation --------------------------------------
anno <- AnnotationDbi::select(
  hgu133plus2.db,
  keys      = rownames(expr.mat),
  columns   = c("PROBEID", "SYMBOL"),
  keytype   = "PROBEID"
)
anno <- anno[!is.na(anno$SYMBOL) & anno$SYMBOL != "", ]
# select() may return a 1:many mapping; keep the first symbol per probe.
anno <- anno[!duplicated(anno$PROBEID), , drop = FALSE]
expr.mat <- expr.mat[anno$PROBEID, , drop = FALSE]
sym <- setNames(anno$SYMBOL, anno$PROBEID)

# --- 2) Collapse multiple probes -> median per gene --------------------------
# Group probes BY GENE SYMBOL (not by expression value) and take the median
# across probes within each group, separately for every sample.
grp  <- split(seq_along(sym), unname(sym))
vals <- lapply(grp, function(ii) {
  m <- expr.mat[ii, , drop = FALSE]
  if (nrow(m) == 1L) as.numeric(m[1, ]) else apply(m, 2, median, na.rm = TRUE)
})
proc <- matrix(unlist(vals, use.names = FALSE), nrow = length(grp), byrow = TRUE,
               dimnames = list(names(grp), colnames(expr.mat)))
proc <- proc[order(rownames(proc)), , drop = FALSE]
cat(sprintf("[01] Annotated + collapsed to %d unique genes x %d samples.\n",
            nrow(proc), ncol(proc)))

# --- 3) Identify lesional/non-lesional status and patient pairing ------------
# GSE30999 titles are "<subject-id>_NL" / "<subject-id>_LS"; the
# "biopsy type:ch1" and "subject:ch1" metadata columns are used to cross-check
# this parsing rather than relying on a free-text keyword match.
ttl <- as.character(pdat[["title"]])
typ <- ifelse(grepl("_NL$", ttl), "NL",
       ifelse(grepl("_LS$", ttl), "LS", NA_character_))
pid <- sub("_(NL|LS)$", "", ttl)
sample_id <- rownames(pdat)

if ("biopsy type:ch1" %in% colnames(pdat)) {
  bt <- as.character(pdat[["biopsy type:ch1"]])
  bad <- !is.na(typ) & ((typ == "NL" & bt != "non-lesion (NL)") |
                        (typ == "LS" & bt != "psoriasis lesion (LS)"))
  if (any(bad)) stop("[01] Lesion label in title conflicts with biopsy type:ch1.")
  cat("[01] Lesion labels cross-checked against 'biopsy type:ch1'.\n")
}
if ("subject:ch1" %in% colnames(pdat)) {
  s2 <- as.character(pdat[["subject:ch1"]])
  if (!identical(s2, pid))
    warning("[01] 'subject:ch1' differs from the title prefix; using title prefix.")
  else cat("[01] Patient IDs cross-checked against 'subject:ch1'.\n")
}
if (anyNA(typ)) stop("[01] Unparseable sample title(s); inspect pdat$title.")

les <- split(sample_id[typ == "LS"], pid[typ == "LS"])
non <- split(sample_id[typ == "NL"], pid[typ == "NL"])
common <- sort(intersect(names(les), names(non)))
pairs_df <- data.frame(
  patient_id         = common,
  lesional_sample    = unlist(les[common], use.names = FALSE),
  nonlesional_sample = unlist(non[common], use.names = FALSE),
  stringsAsFactors   = FALSE
)
rownames(pairs_df) <- NULL
cat(sprintf("[01] %d subjects with an NL sample, %d with an LS sample; %d complete pairs.\n",
            length(non), length(les), nrow(pairs_df)))

cat(sprintf("[01] Identified %d complete lesional/non-lesional pairs (manuscript target = 81).\n",
            nrow(pairs_df)))
if (nrow(pairs_df) != 81)
  warning("[01] Pair count != 81; re-check the lesion-status / patient-ID parsing above.")

saveRDS(proc,      file.path("data", "gse30999_processed.rds"))
saveRDS(pairs_df,  file.path("data", "gse30999_pairs.rds"))
cat("[01] Saved processed matrix + pairing table to data/.\n")
