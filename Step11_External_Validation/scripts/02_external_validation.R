#!/usr/bin/env Rscript
# =============================================================================
# 02_external_validation.R
# -----------------------------------------------------------------------------
# Purpose : External validation of the PPAR-responsive gene signature and the
#           four hub genes in the independent cohort GSE30999.
#
# Methods (aligned with manuscript Methods 2.2.11):
#   - Differential expression between PAIRED lesional and non-lesional skin:
#       Wilcoxon signed-rank test, Benjamini-Hochberg (BH) FDR applied WITHIN
#       the pre-specified 18-gene set (family-wise correction inside the set,
#       not genome-wide).
#   - Discriminative performance: ROC / AUC, 95% CI by the DeLong method.
#   - Direction concordance: GSE30999 log2FC sign vs discovery-cohort direction.
#
# Input  : data/gse30999_processed.rds, data/gse30999_pairs.rds
#          scripts/input/signature_18_genes.csv
#          scripts/input/hub_genes.csv
# Output : results/external_validation_table.csv (same columns as the
#          supplementary external-validation table shipped with the manuscript)
#          figures/Fig16_external_validation.pdf     (A: paired boxplots,
#                                                      B: ROC curves of 4 hubs)
#          figures/Fig16_external_validation.png (300 dpi; see below for TIFF)
#
# Usage  : Rscript scripts/02_external_validation.R
#
# NOTE: requires R >= 4.3 with pROC, ggplot2, reshape2 and gridExtra.
#       ROC direction is resolved automatically (positive class = non-lesional)
#       and AUC confidence intervals use pROC::ci.auc (DeLong).
# =============================================================================

suppressMessages({
  if (!requireNamespace("pROC", quietly = TRUE))      install.packages("pROC")
  if (!requireNamespace("ggplot2", quietly = TRUE))   install.packages("ggplot2")
  if (!requireNamespace("reshape2", quietly = TRUE))  install.packages("reshape2")
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  library(pROC); library(ggplot2); library(reshape2); library(gridExtra)
})

proc     <- readRDS("data/gse30999_processed.rds")
pairs_df <- readRDS("data/gse30999_pairs.rds")
sig      <- read.csv("scripts/input/signature_18_genes.csv", stringsAsFactors = FALSE)
hub      <- read.csv("scripts/input/hub_genes.csv",          stringsAsFactors = FALSE)

# --- Subset expression to the 81 paired samples, ordered by pair -------------
genes_in <- intersect(sig$Gene, rownames(proc))
missing  <- setdiff(sig$Gene, rownames(proc))
if (length(missing) > 0)
  warning(sprintf("[02] %d signature gene(s) absent after probe annotation: %s",
                  length(missing), paste(missing, collapse = ", ")))

les_mat <- proc[genes_in, pairs_df$lesional_sample,    drop = FALSE]
non_mat <- proc[genes_in, pairs_df$nonlesional_sample, drop = FALSE]
n_pairs <- nrow(pairs_df)

# --- Per-gene validation -----------------------------------------------------
rows <- list()
roc_list <- list()
p_raw <- numeric(0)
for (g in genes_in) {
  x <- as.numeric(les_mat[g, ])   # lesional
  y <- as.numeric(non_mat[g, ])   # non-lesional
  wt <- wilcox.test(x, y, paired = TRUE)
  p  <- wt$p.value
  p_raw[g] <- p

  med_l  <- median(x); med_n <- median(y)
  log2fc <- med_l - med_n

  scores <- c(x, y)
  resp   <- c(rep(0, length(x)), rep(1, length(y)))   # 1 = non-lesional (positive class)
  roc_obj <- roc(response = resp, predictor = scores,
                 direction = "auto", quiet = TRUE)
  auc_raw <- as.numeric(roc_obj$auc)
  ci_obj  <- as.numeric(ci.auc(roc_obj, method = "delong", quiet = TRUE))
  # Report the oriented discrimination AUC (max[AUC, 1-AUC]) so that the AUC
  # and its DeLong CI always refer to the same direction as the column name
  # "non-lesional vs lesional". direction="auto" can return AUC < 0.5 when
  # ties are present; flip the CI endpoints together with the AUC.
  if (auc_raw < 0.5) {
    auc_val <- 1 - auc_raw
    ci_obj  <- c(1 - ci_obj[3], 1 - ci_obj[2], 1 - ci_obj[1])
  } else {
    auc_val <- auc_raw
  }
  disc    <- auc_val

  disc_dir <- sig$discovery_direction[sig$Gene == g]
  consistent <- ifelse(disc_dir == "up", log2fc > 0,
                ifelse(disc_dir == "down", log2fc < 0, NA))

  rows[[g]] <- data.frame(
    Gene = g,
    Direction_in_discovery_cohort = disc_dir,
    Median_expression_lesional   = round(med_l, 3),
    Median_expression_nonlesional = round(med_n, 3),
    log2FC_lesional_over_nonlesional = round(log2fc, 3),
    AUC_nonlesional_vs_lesional  = round(auc_val, 4),
    AUC_95CI_DeLong = sprintf("%.3f-%.3f", ci_obj[1], ci_obj[3]),
    Discrimination_maxAUC_1minusAUC = round(disc, 4),
    P_value_Wilcoxon_paired = format(p, scientific = TRUE, digits = 3),
    Direction_consistent = consistent,
    stringsAsFactors = FALSE
  )

  roc_list[[g]] <- data.frame(
    Gene = g,
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    AUC = round(auc_val, 3),
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, rows)
res$BH_adjusted_FDR <- p.adjust(p_raw[rownames(res)], method = "BH")
res$BH_adjusted_FDR <- format(res$BH_adjusted_FDR, scientific = TRUE, digits = 3)
p_raw <- p_raw[order(names(p_raw))]
res <- res[order(res$Gene), ]
rownames(res) <- NULL

# --- Save S4-style table -----------------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
out_csv <- "results/external_validation_table.csv"
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("[02] Wrote %s (%d genes).\n", out_csv, nrow(res)))

# --- Summary (should match manuscript: 17/18 direction concordant,
#     18/18 FDR < 0.05) -------------------------------------------------------
fdr_num <- as.numeric(res$BH_adjusted_FDR)
n_sig   <- sum(fdr_num < 0.05, na.rm = TRUE)
n_cons  <- sum(res$Direction_consistent == TRUE, na.rm = TRUE)
cat(sprintf("[02] Genes with BH-FDR < 0.05: %d / %d\n", n_sig, nrow(res)))
cat(sprintf("[02] Directionally concordant:   %d / %d\n", n_cons, nrow(res)))

# =============================================================================
# Figure 16: external validation of the FOUR hub genes in GSE30999
#   (A) paired boxplots: non-lesional (control) vs lesional (disease)
#   (B) ROC curves with AUC (95% CI) labels
# =============================================================================
hub_genes <- intersect(hub$Gene, genes_in)
if (length(hub_genes) != 4)
  warning(sprintf("[02] Expected 4 hub genes; found %d: %s",
                  length(hub_genes), paste(hub_genes, collapse = ", ")))

# --- Fig 16A: paired boxplots -----------------------------------------------
# melt() on a matrix returns Var1 (row = pair), Var2 (col = gene), value.
les_long <- reshape2::melt(t(les_mat[hub_genes, , drop = FALSE]))
non_long <- reshape2::melt(t(non_mat[hub_genes, , drop = FALSE]))
box_df <- data.frame(Pair        = les_long$Var1,
                     Gene        = les_long$Var2,
                     Lesional    = les_long$value,
                     Nonlesional = non_long$value)
box_long <- reshape2::melt(box_df, id.vars = c("Gene", "Pair"),
                           variable.name = "Group", value.name = "Expression")
box_long$Group <- ifelse(box_long$Group == "Lesional", "Lesional", "Non-lesional")

pA <- ggplot(box_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Non-lesional" = "#2C7BB6", "Lesional" = "#D7191C")) +
  labs(x = NULL, y = "Expression (RMA, log2)",
       title = "Fig 16A. Hub-gene expression in paired skin (GSE30999, n = 81 pairs)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey90"))

# --- Fig 16B: ROC curves -----------------------------------------------------
roc_all <- do.call(rbind, roc_list[hub_genes])
lab_tab <- res[res$Gene %in% hub_genes, c("Gene", "AUC_nonlesional_vs_lesional",
                                          "AUC_95CI_DeLong")]
lab_tab$lab <- sprintf("%s (AUC %s, 95%%CI %s)", lab_tab$Gene,
                       lab_tab$AUC_nonlesional_vs_lesional, lab_tab$AUC_95CI_DeLong)

pB <- ggplot(roc_all, aes(x = FPR, y = TPR, colour = Gene)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_brewer(palette = "Set1",
                      labels = setNames(lab_tab$lab, lab_tab$Gene)) +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = "Hub gene",
       title = "Fig 16B. ROC of four hub genes (non-lesional vs lesional)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

fig_pdf <- "figures/Fig16_external_validation.pdf"
fig_tif <- "figures/Fig16_external_validation.tif"
pdf(fig_pdf, width = 11, height = 4.5)
gridExtra::grid.arrange(pA, pB, nrow = 1)
dev.off()

# Portable raster output: the default macOS quartz TIFF device ignores
# compression="lzw" (producing a ~17 MB file) and the cairo device requires
# X11 libraries that are not present on every machine. PDF + 300-dpi PNG are
# therefore written here; a journal-ready LZW TIFF can be obtained with:
#   python3 -c "from PIL import Image; im = Image.open('figures/Fig16_external_validation.png'); \\
#     im.save('figures/Fig16_external_validation.tif', compression='tiff_lzw', dpi=(300,300))"
# A 300-dpi LZW TIFF reference copy is shipped in results/.
fig_png <- "figures/Fig16_external_validation.png"
png(fig_png, width = 28, height = 11.5, units = "cm", res = 300)
gridExtra::grid.arrange(pA, pB, nrow = 1)
dev.off()
cat(sprintf("[02] Wrote %s and %s.\n", fig_pdf, fig_png))
