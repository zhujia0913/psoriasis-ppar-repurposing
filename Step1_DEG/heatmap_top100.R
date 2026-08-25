#!/usr/bin/env Rscript
# Heatmap: Top 50 up + Top 50 down DEGs (Z-score)
# Publication-quality output for PLoS One
# Author: auto-generated 2026-06-19

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  script_path <- sub("^--file=", "", file_arg)
} else {
  here <- getwd()
  script_path <- "heatmap_top100.R"
}

# ── Paths ──────────────────────────────────────────────────────────────
expr_file <- file.path(here, "expression_combat_corrected.csv")
deg_file  <- file.path(here, "deg_results_significant.tsv")
samp_file <- file.path(here, "sample_info.tsv")
out_dir   <- file.path(here, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Read data ───────────────────────────────────────────────────────
message("Reading expression matrix...")
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)

message("Reading DEG table...")
deg <- read.delim(deg_file, stringsAsFactors = FALSE)

message("Reading sample info...")
si <- read.delim(samp_file, stringsAsFactors = FALSE)
rownames(si) <- si$sample
si <- si[colnames(expr), ]

# ── 2. Select top 50 up / top 50 down ──────────────────────────────────
deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
deg_sig <- deg_sig[order(deg_sig$logFC, decreasing = TRUE), ]

top_up   <- head(deg_sig, 50)
top_down <- tail(deg_sig, 50)
top_all  <- rbind(top_up, top_down)

genes_use <- top_all$Gene
message(sprintf("Selected %d genes (50 up + 50 down)", length(genes_use)))

# ── 3. Subset expression + Z-score ─────────────────────────────────────
mat <- as.matrix(expr[genes_use, , drop = FALSE])
mat_z <- t(scale(t(mat)))                     # row-wise Z-score
mat_z <- pmin(pmax(mat_z, -3), 3)

# ── 4. Column ordering: Group(NN→PP) → Dataset ───────────────────────
si$dataset <- si$batch
ord <- order(factor(si$group,   levels = c("NN", "PP")),
             factor(si$dataset, levels = c("GSE13355", "GSE14905", "GSE78097")))
mat_z  <- mat_z[, ord, drop = FALSE]
si_ord <- si[ord, , drop = FALSE]

# ── 5. Annotation bars (Group on top, Dataset below) ──────────────────
ann_col <- data.frame(
  Group   = si_ord$group,
  Dataset = si_ord$dataset,
  row.names = rownames(si_ord)
)
ann_colors <- list(
  Group   = c(NN = "#FDBF6F", PP = "#A6CEE3"),
  Dataset = c(GSE13355 = "#E41A1C", GSE14905 = "#377EB8", GSE78097 = "#4DAF4A")
)

# ── 6. Color palette (blue=low, white=mid, red=high) ──────────────────
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ── 7. Column gaps: only group-level (NN→PP) ─────────────────────────
nn_cnt <- sum(si_ord$group == "NN")
gaps_col <- nn_cnt  # single gap between NN and PP

# ── 9. Heatmap ─────────────────────────────────────────────────────────
# Row font size: auto-scale based on gene count
n_genes <- length(genes_use)
fs_row  <- if (n_genes <= 50) 8 else if (n_genes <= 100) 7 else 6

png(file.path(out_dir, "Heatmap_Top100_DEG.png"),
    width = 5400, height = 4800, res = 300)

pheatmap(mat_z,
  color             = heat_colors,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = fs_row,
  fontsize          = 11,
  gaps_col          = gaps_col,
  main              = "Top 50 Up-regulated and 50 Down-regulated DEGs",
  treeheight_row    = 60,
  legend_breaks     = seq(-3, 3, 1),
  legend_labels     = seq(-3, 3, 1),
  border_color      = NA
)
dev.off()
message("PNG saved.")

pdf(file.path(out_dir, "Heatmap_Top100_DEG.pdf"), width = 18, height = 16)
pheatmap(mat_z,
  color             = heat_colors,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 6,
  fontsize          = 9,
  gaps_col          = gaps_col,
  main              = "Top 50 Up-regulated and 50 Down-regulated DEGs",
  treeheight_row    = 60,
  legend_breaks     = seq(-3, 3, 1),
  legend_labels     = seq(-3, 3, 1),
  border_color      = NA
)
dev.off()
message("PDF saved.")

# ── 10. Source & gene list ─────────────────────────────────────────────
writeLines(
  c(sprintf("Script: %s", script_path),
    sprintf("Expression: expression_combat_corrected.csv (%d x %d)",
            nrow(expr), ncol(expr)),
    "DEG: deg_results_significant.tsv",
    "Sample info: sample_info.tsv",
    sprintf("Genes: %d (top 50 by logFC + bottom 50 by logFC)", n_genes),
    "Z-score capped at ±3, RdBu divergent palette (blue=low, red=high)",
    "Columns ordered: NN(GSE13355→GSE14905→GSE78097) → PP(GSE13355→GSE14905→GSE78097)",
    "Annotation: Group(top) + Dataset(bottom), row Direction bar removed",
    "Gaps: group-level (NN→PP) + dataset-level within each group"),
  file.path(out_dir, "Heatmap_Top100_DEG.source")
)

# Save gene list too
write.table(top_all, file.path(out_dir, "top100_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Done.")
