#!/usr/bin/env Rscript
# =============================================================================
# 02_external_validation.R
# -----------------------------------------------------------------------------
# Purpose : External validation of the PPAR-responsive 18-gene signature and
#           the four hub genes in the independent cohort GSE30999, and
#           regeneration of the published table and Figure 16.
#
# Methods (manuscript Methods 2.2.11):
#   - Input: ArrayExpress E-GEOD-30999 processed (RMA, log2) values.
#   - Probe collapse: MaxMean (one representative probe per gene, the probe
#     with the highest mean expression across all 170 samples), identical to
#     the discovery-cohort preprocessing in Step 1. The selected probes for
#     the 18 reported genes are shipped in
#     scripts/input/representative_probes.csv and in the supplementary table.
#   - Differential expression between PAIRED lesional and non-lesional skin:
#     Wilcoxon signed-rank test; Benjamini-Hochberg FDR applied WITHIN the
#     pre-specified 18-gene set.
#   - Discriminative performance: ROC / AUC, 95% CI by the DeLong method;
#     positive class = non-lesional, direction resolved automatically.
#   - Direction concordance: GSE30999 log2FC sign vs discovery-cohort direction.
#
# Input  : data/gse30999_processed.rds, data/gse30999_pairs.rds (from 01)
#          scripts/input/signature_18_genes.csv, hub_genes.csv
# Output : results/external_validation_table.csv   (same numbers as the xlsx)
#          figures/Fig16_external_validation.pdf / .png (600 dpi, 5908x5413)
#
# Usage  : Rscript scripts/02_external_validation.R
# Requires: R >= 4.3, pROC, ggplot2, gridExtra, ragg (all CRAN).
# =============================================================================

suppressMessages({
  for (p in c("pROC", "ggplot2", "gridExtra", "ragg"))
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(pROC); library(ggplot2); library(gridExtra)
})

proc     <- readRDS("data/gse30999_processed.rds")
pairs_df <- readRDS("data/gse30999_pairs.rds")
sig      <- read.csv("scripts/input/signature_18_genes.csv", stringsAsFactors = FALSE)
hubdf    <- read.csv("scripts/input/hub_genes.csv", stringsAsFactors = FALSE)
sel      <- read.csv("scripts/input/representative_probes.csv", stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

genes_in <- intersect(sig$Gene, rownames(proc))
if (length(genes_in) != 18L)
  stop(sprintf("[02] Expected 18 signature genes, found %d", length(genes_in)))

les_mat <- proc[genes_in, pairs_df$lesional_sample,    drop = FALSE]
non_mat <- proc[genes_in, pairs_df$nonlesional_sample, drop = FALSE]
cat(sprintf("[02] %d complete pairs; %d genes\n", nrow(pairs_df), length(genes_in)))

# --- Per-gene statistics -----------------------------------------------------
rows <- list(); p_raw <- numeric(0); roc_list <- list()
for (g in genes_in) {
  x <- as.numeric(les_mat[g, ]); y <- as.numeric(non_mat[g, ])
  p <- wilcox.test(x, y, paired = TRUE)$p.value
  p_raw[g] <- p
  med_l <- median(x); med_n <- median(y); log2fc <- med_l - med_n

  scores <- c(x, y)
  resp   <- c(rep(0, length(x)), rep(1, length(y)))
  roc_obj <- roc(response = resp, predictor = scores,
                 direction = "auto", quiet = TRUE)
  auc_raw <- as.numeric(roc_obj$auc)
  ci_obj  <- as.numeric(ci.auc(roc_obj, method = "delong", quiet = TRUE))
  # Keep AUC and CI on the same ROC direction (direction="auto" can return
  # AUC < 0.5 when ties are present).
  if (auc_raw < 0.5) {
    auc_val <- 1 - auc_raw
    ci_obj  <- c(1 - ci_obj[3], 1 - ci_obj[2], 1 - ci_obj[1])
  } else auc_val <- auc_raw

  disc_dir   <- sig$discovery_direction[sig$Gene == g]
  consistent <- ifelse(disc_dir == "up", log2fc > 0,
                ifelse(disc_dir == "down", log2fc < 0, NA))

  rows[[g]] <- data.frame(
    Gene = g,
    Probe_ID_MaxMean = sel$Probe[sel$Gene == g],
    Direction_in_discovery_cohort = disc_dir,
    Median_expression_lesional = med_l,
    Median_expression_nonlesional = med_n,
    log2FC_lesional_over_nonlesional = log2fc,
    AUC_nonlesional_vs_lesional = auc_val,
    AUC_95CI_DeLong = sprintf("%.3f\u2013%.3f", ci_obj[1], ci_obj[3]),
    Discrimination_maxAUC_1minusAUC = auc_val,
    P_value_Wilcoxon_paired = p,
    Direction_consistent = consistent,
    stringsAsFactors = FALSE)
  roc_list[[g]] <- roc_obj
}
res <- do.call(rbind, rows)
res$BH_adjusted_FDR <- p.adjust(p_raw[res$Gene], method = "BH")
res <- res[, c("Gene", "Probe_ID_MaxMean", "Direction_in_discovery_cohort",
               "Median_expression_lesional", "Median_expression_nonlesional",
               "log2FC_lesional_over_nonlesional", "AUC_nonlesional_vs_lesional",
               "AUC_95CI_DeLong", "Discrimination_maxAUC_1minusAUC",
               "P_value_Wilcoxon_paired", "BH_adjusted_FDR",
               "Direction_consistent")]
res <- res[order(match(res$Gene, sig$Gene)), ]; rownames(res) <- NULL
write.csv(res, "results/external_validation_table.csv", row.names = FALSE)
cat(sprintf("[02] Wrote results/external_validation_table.csv (%d genes)\n", nrow(res)))
cat(sprintf("[02] BH-FDR < 0.05: %d/18 | directionally concordant: %d/18\n",
            sum(res$BH_adjusted_FDR < 0.05),
            sum(res$Direction_consistent == TRUE, na.rm = TRUE)))

# =============================================================================
# Figure 16 (journal layout): A = paired boxplots, B = ROC of four hub genes
# =============================================================================
# geom_jitter() consumes random numbers at draw time; fix the seed so the
# archived TIFF can be regenerated pixel-for-pixel.
set.seed(20261111)
hub <- intersect(hubdf$Gene, genes_in)
hub <- c("ADIPOQ", "APOE", "PLIN1", "BCL2")[c("ADIPOQ", "APOE", "PLIN1", "BCL2") %in% genes_in]
box_df <- do.call(rbind, lapply(hub, function(g) {
  rbind(data.frame(Gene = g, Group = "control",
                   Expression = as.numeric(non_mat[g, ])),
        data.frame(Gene = g, Group = "disease",
                   Expression = as.numeric(les_mat[g, ])))
}))
box_df$Group <- factor(box_df$Group, levels = c("control", "disease"))
box_df$Gene  <- factor(box_df$Gene,  levels = hub)

roc_df <- do.call(rbind, lapply(hub, function(g) {
  ro <- roc_list[[g]]
  data.frame(Gene = g, FPR = 1 - ro$specificities, TPR = ro$sensitivities)
}))
lab <- do.call(rbind, lapply(hub, function(g) {
  r <- res[res$Gene == g, ]
  ci <- strsplit(r$AUC_95CI_DeLong, "\u2013")[[1]]
  data.frame(Gene = g,
             label = sprintf("%s (AUC = %.3f, 95%% CI %s\u2013%s)",
                             g, r$AUC_nonlesional_vs_lesional, ci[1], ci[2]))
}))

col_ctrl <- "#2C7BB6"; col_dis <- "#D7191C"
roc_cols <- c(ADIPOQ = "#D7191C", APOE = "#8E44AD",
              PLIN1 = "#1F3B8B", BCL2 = "#2CA02C")

pA <- ggplot(box_df, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, size = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.22, dodge.width = 0.62),
              size = 0.6, alpha = 0.5) +
  scale_fill_manual(values = c(control = col_ctrl, disease = col_dis),
                    labels = c("control", "disease")) +
  labs(title = "Expression of Hub Genes in GSE30999",
       subtitle = "NL: Control, LS: Disease",
       x = "Genes", y = "Normalized Expression Level (log2)",
       fill = "Group", tag = "A") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0, size = 11.5),
        plot.tag = element_text(face = "bold", size = 20),
        legend.position = "top", legend.title = element_text(size = 11.5),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12.5), axis.text = element_text(size = 11.5),
        plot.margin = unit(c(8, 10, 8, 20), "pt"),
        panel.grid = element_line(linewidth = 0.4))
ann <- data.frame(Gene = factor(hub, levels = hub),
                  y = sapply(hub, function(g) max(box_df$Expression[box_df$Gene == g]) * 1.06),
                  lab = rep("p < 0.001", length(hub)))
pA <- pA + geom_text(data = ann, aes(x = Gene, y = y, label = lab),
                     inherit.aes = FALSE, size = 3.8)

pB <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Gene)) +
  geom_line(linewidth = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey55", linewidth = 0.6) +
  scale_colour_manual(values = roc_cols, breaks = hub,
                      labels = setNames(lab$label, lab$Gene)) +
  labs(title = "ROC Curves of Hub Genes",
       subtitle = "Diagnostic Performance in GSE30999 Dataset",
       x = "1 \u2013 Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)", colour = "Genes", tag = "B") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 11.5),
        plot.tag = element_text(face = "bold", size = 20),
        legend.position = "right", legend.title = element_text(size = 11.5),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12.5), axis.text = element_text(size = 11.5),
        plot.margin = unit(c(8, 10, 8, 20), "pt"),
        panel.grid = element_line(linewidth = 0.4))

pdf("figures/Fig16_external_validation.pdf", width = 9.85, height = 9.02)
grid.arrange(pA, pB, nrow = 2, heights = c(1, 1)); dev.off()
ragg::agg_png("figures/Fig16_external_validation.png",
              width = 5908, height = 5413, units = "px", res = 600,
              background = "white")
grid.arrange(pA, pB, nrow = 2, heights = c(1, 1)); dev.off()
cat("[02] Wrote figures/Fig16_external_validation.pdf and .png (600 dpi).\n")
cat("[02] For a journal-ready LZW TIFF see figures/README.md.\n")
