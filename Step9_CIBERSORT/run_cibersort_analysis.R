# =============================================================================
# CIBERSORT 免疫浸润分析
# 对齐 Liang et al. 2023 (J Orthop Surg Res) 方法学
# 输入: Pre-ComBat log2表达矩阵 (Affymetrix HG-U133 Plus 2.0 RMA标准化)
# 参数: perm=1000, QN=TRUE (芯片数据需分位数归一化)
# 输出: 22种免疫细胞比例矩阵 + 4张图 + 5个TSV
# =============================================================================

library(CIBERSORT)
cat(sprintf("CIBERSORT包版本: %s\n", as.character(packageVersion("CIBERSORT"))))
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# --- 1. 路径设置 ---
base_dir <- here
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
scripts_dir <- file.path(base_dir, "scripts")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# --- 2. 读取数据 ---
# Pre-ComBat表达矩阵（CIBERSORT要求：基因在行，样本在列，第一列为基因名）
expr_pre <- read.csv(file.path(data_dir, "expression_pre_combat.csv"),
                     row.names = 1, check.names = FALSE)
cat(sprintf("Pre-ComBat表达矩阵: %d基因 x %d样本\n", nrow(expr_pre), ncol(expr_pre)))

# 样本信息
sample_info <- read.delim(file.path(data_dir, "sample_info.tsv"))
names(sample_info)[1:2] <- c("Sample", "Group")
cat(sprintf("样本信息: %d样本 (NN=%d, PP=%d)\n",
            nrow(sample_info),
            sum(sample_info$Group == "NN"),
            sum(sample_info$Group == "PP")))

# 取共有样本
common_samples <- intersect(colnames(expr_pre), sample_info$Sample)
cat(sprintf("共有样本: %d\n", length(common_samples)))

expr_pre <- expr_pre[, common_samples]
sample_info <- sample_info[match(common_samples, sample_info$Sample), ]

# --- 3. 准备CIBERSORT输入文件 ---
# CIBERSORT包的cibersort()函数接受文件路径，不是矩阵对象
# 需要写入临时文件：第一列为基因名，其余为样本表达值
input_file <- file.path(data_dir, "cibersort_input_pre_combat.txt")
expr_df <- data.frame(GeneSymbol = rownames(expr_pre), expr_pre, check.names = FALSE)
write.table(expr_df, file = input_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("CIBERSORT输入文件: %s (%d基因 x %d样本)\n", input_file, nrow(expr_df), ncol(expr_df) - 1))

# LM22签名矩阵
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
cat(sprintf("LM22签名: %s\n", sig_matrix))

# --- 4. 运行CIBERSORT ---
# perm=1000: 置换检验1000次（对齐Liang 2023）
# QN=TRUE: 分位数归一化（Affymetrix芯片数据）
cat("开始运行CIBERSORT (perm=1000, QN=TRUE)...\n")
cat("预计耗时: 5-15分钟\n")

set.seed(42)
cibersort_results <- cibersort(sig_matrix, input_file, perm = 1000, QN = TRUE)

cat(sprintf("CIBERSORT完成: %d样本 x %d列\n", nrow(cibersort_results), ncol(cibersort_results)))

# --- 5. 保存原始结果 ---
# 结果包含22种细胞比例 + P-value + Correlation + RMSE
output_raw <- file.path(results_dir, "cibersort_results_raw.tsv")
write.table(cibersort_results, file = output_raw, sep = "\t", quote = FALSE)
cat(sprintf("原始结果保存: %s\n", output_raw))

# --- 6. 过滤显著样本 (P < 0.05, 对齐Liang 2023) ---
significant_samples <- rownames(cibersort_results)[cibersort_results[, "P-value"] < 0.05]
cat(sprintf("显著样本 (P<0.05): %d / %d (%.1f%%)\n",
            length(significant_samples), nrow(cibersort_results),
            100 * length(significant_samples) / nrow(cibersort_results)))

# 提取22种细胞比例（只保留显著样本）
cell_cols <- colnames(cibersort_results)[1:22]
results_filtered <- cibersort_results[significant_samples, cell_cols, drop = FALSE]
results_all <- cibersort_results[, cell_cols, drop = FALSE]

# 保存过滤后和完整结果
write.table(results_filtered, file.path(results_dir, "cibersort_proportions_filtered.tsv"),
            sep = "\t", quote = FALSE)
write.table(results_all, file.path(results_dir, "cibersort_proportions_all.tsv"),
            sep = "\t", quote = FALSE)

# --- 7. 合并分组信息 ---
results_with_group <- data.frame(
  Sample = rownames(results_all),
  Group = sample_info$Group[match(rownames(results_all), sample_info$Sample)],
  results_all,
  `P-value` = cibersort_results[, "P-value"],
  Correlation = cibersort_results[, "Correlation"],
  RMSE = cibersort_results[, "RMSE"],
  check.names = FALSE
)
write.table(results_with_group, file.path(results_dir, "cibersort_results_with_group.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 8. 图1: 免疫细胞浸润热图 (22种细胞 × 样本) ---
# 对齐Liang 2023 Fig 6A
annotation_col <- data.frame(
  Group = factor(sample_info$Group[match(rownames(results_all), sample_info$Sample)],
                 levels = c("NN", "PP")),
  row.names = rownames(results_all)
)
ann_colors <- list(Group = c(NN = "#4DBBD5", PP = "#E64B35"))

png(file.path(figures_dir, "heatmap_cibersort_infiltration.png"),
    width = 2400, height = 800, res = 150)
pheatmap(
  t(results_all),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  fontsize_row = 9,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "CIBERSORT Immune Cell Infiltration (22 cell types)"
)
dev.off()
cat("图1保存: heatmap_cibersort_infiltration.png\n")

# --- 9. 图2: 22种免疫细胞间Spearman相关性热图 ---
# 对齐Liang 2023 Fig 6B
cor_matrix <- cor(results_all, method = "spearman")

# 自定义计算p值矩阵
cor_pmat_custom <- function(mat, method = "spearman") {
  n <- ncol(mat)
  pmat <- matrix(NA, n, n, dimnames = list(colnames(mat), colnames(mat)))
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        pmat[i, j] <- 0
      } else {
        pmat[i, j] <- cor.test(mat[, i], mat[, j], method = method)$p.value
      }
    }
  }
  return(pmat)
}

png(file.path(figures_dir, "heatmap_cibersort_correlation.png"),
    width = 1400, height = 1200, res = 150)
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.8,
  col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(100),
  title = "Spearman Correlation among 22 Immune Cell Types",
  mar = c(0, 0, 2, 0),
  addCoef.col = "black",
  number.cex = 0.5,
  p.mat = cor_pmat_custom(results_all),
  sig.level = 0.05,
  insig = "blank"
)
dev.off()
cat("图2保存: heatmap_cibersort_correlation.png\n")

# --- 10. 图3: 组间差异小提琴图 ---
# 对齐Liang 2023 Fig 6C (violin plot)
results_long <- results_with_group %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "CellType",
    values_to = "Proportion"
  )

# Wilcoxon rank-sum test
pvals <- sapply(cell_cols, function(ct) {
  nn_vals <- results_all[results_with_group$Group == "NN", ct]
  pp_vals <- results_all[results_with_group$Group == "PP", ct]
  wilcox.test(nn_vals, pp_vals)$p.value
})
padj <- p.adjust(pvals, method = "BH")
sig_labels <- ifelse(padj < 0.001, "***",
                     ifelse(padj < 0.01, "**",
                            ifelse(padj < 0.05, "*", "ns")))

diff_df <- data.frame(
  CellType = cell_cols,
  P.value = pvals,
  P.adj = padj,
  Significance = sig_labels,
  stringsAsFactors = FALSE
)
diff_df <- diff_df[order(diff_df$P.adj), ]
write.table(diff_df, file.path(results_dir, "cibersort_comparison_NN_vs_PP.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 小提琴图
results_long$CellType <- factor(results_long$CellType, levels = cell_cols)
results_long$Group <- factor(results_long$Group, levels = c("NN", "PP"))

png(file.path(figures_dir, "violin_cibersort_comparison.png"),
    width = 2400, height = 1200, res = 150)
p <- ggviolin(
  results_long,
  x = "Group",
  y = "Proportion",
  fill = "Group",
  palette = c("#4DBBD5", "#E64B35"),
  add = "boxplot",
  add.params = list(width = 0.1, fill = "white"),
  facet.by = "CellType",
  scales = "free_y",
  short.panel.labs = TRUE,
  panel.labs.background = list(fill = "white", color = "grey80")
) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  labs(title = "Immune Cell Infiltration: NN vs PP (CIBERSORT)",
       x = NULL, y = "Proportion") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 7),
    legend.position = "bottom"
  )
print(p)
dev.off()
cat("图3保存: violin_cibersort_comparison.png\n")

# --- 11. 图4: Hub基因 × 免疫细胞Spearman相关性热图 ---
# 我们独有的增量分析（Liang 2023未做）
hub_genes <- c("ADIPOQ", "APOE", "PLIN1", "BCL2")
hub_expr <- expr_pre[hub_genes, , drop = FALSE]

# 计算Hub基因与22种免疫细胞的Spearman相关
cor_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols),
                  dimnames = list(hub_genes, cell_cols))
p_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols),
                dimnames = list(hub_genes, cell_cols))

for (gene in hub_genes) {
  for (ct in cell_cols) {
    gene_vec <- as.numeric(hub_expr[gene, ])
    test <- cor.test(gene_vec, results_all[, ct], method = "spearman")
    cor_hub[gene, ct] <- test$estimate
    p_hub[gene, ct] <- test$p.value
  }
}

# FDR校正
p_hub_adj <- matrix(p.adjust(as.vector(p_hub), method = "BH"),
                    nrow = length(hub_genes), ncol = length(cell_cols),
                    dimnames = list(hub_genes, cell_cols))

# 保存相关性结果
hub_cor_df <- data.frame(
  Gene = rep(hub_genes, each = length(cell_cols)),
  CellType = rep(cell_cols, length(hub_genes)),
  Rho = as.vector(cor_hub),
  P.value = as.vector(p_hub),
  P.adj = as.vector(p_hub_adj)
)
hub_cor_df <- hub_cor_df[order(hub_cor_df$P.adj), ]
write.table(hub_cor_df, file.path(results_dir, "hub_immune_correlation_cibersort.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 热图
sig_stars <- ifelse(p_hub_adj < 0.001, "***",
                    ifelse(p_hub_adj < 0.01, "**",
                           ifelse(p_hub_adj < 0.05, "*", "")))

png(file.path(figures_dir, "heatmap_hub_immune_cibersort.png"),
    width = 1800, height = 800, res = 150)
pheatmap(
  cor_hub,
  display_numbers = sig_stars,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  fontsize_row = 12,
  fontsize_col = 9,
  fontsize_number = 10,
  number_color = "black",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.6, 0.6, length.out = 101),
  main = "Hub Genes × Immune Cells Spearman Correlation (CIBERSORT)"
)
dev.off()
cat("图4保存: heatmap_hub_immune_cibersort.png\n")

# --- 12. 图5: PCA聚类图 (对齐Liang 2023) ---
# 基于免疫细胞浸润矩阵做PCA
pca_result <- prcomp(results_all, scale. = TRUE, center = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = factor(sample_info$Group[match(rownames(results_all), sample_info$Sample)],
                 levels = c("NN", "PP"))
)
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

png(file.path(figures_dir, "pca_cibersort.png"),
    width = 1000, height = 800, res = 150)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.95, type = "t") +
  scale_color_manual(values = c(NN = "#4DBBD5", PP = "#E64B35")) +
  labs(title = "PCA of Immune Cell Infiltration (CIBERSORT)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_bw() +
  theme(legend.position = "bottom")
print(p)
dev.off()
cat("图5保存: pca_cibersort.png\n")

# --- 13. 汇总输出 ---
cat("\n")
cat("========================================\n")
cat("CIBERSORT分析完成\n")
cat("========================================\n")
cat(sprintf("输入: %d基因 x %d样本 (Pre-ComBat)\n", nrow(expr_pre), ncol(expr_pre)))
cat(sprintf("显著样本 (P<0.05): %d / %d\n", length(significant_samples), nrow(cibersort_results)))
cat(sprintf("免疫细胞类型: 22种 (LM22)\n"))
cat(sprintf("\n组间差异 (Wilcoxon, BH-FDR < 0.05):\n"))
sig_cells <- diff_df[diff_df$P.adj < 0.05, ]
if (nrow(sig_cells) > 0) {
  for (i in 1:nrow(sig_cells)) {
    cat(sprintf("  %s: P.adj=%.2e %s\n", sig_cells$CellType[i], sig_cells$P.adj[i], sig_cells$Significance[i]))
  }
} else {
  cat("  无显著差异细胞类型\n")
}
cat(sprintf("\nHub基因显著相关 (P.adj < 0.05): %d / %d\n",
            sum(p_hub_adj < 0.05), length(p_hub_adj)))
cat(sprintf("\n输出文件:\n"))
cat(sprintf("  results/cibersort_results_raw.tsv (原始结果+P/Cor/RMSE)\n"))
cat(sprintf("  results/cibersort_proportions_filtered.tsv (P<0.05过滤)\n"))
cat(sprintf("  results/cibersort_proportions_all.tsv (全部样本22细胞比例)\n"))
cat(sprintf("  results/cibersort_results_with_group.tsv (含分组)\n"))
cat(sprintf("  results/cibersort_comparison_NN_vs_PP.tsv (组间差异统计)\n"))
cat(sprintf("  results/hub_immune_correlation_cibersort.tsv (Hub-免疫相关)\n"))
cat(sprintf("  figures/heatmap_cibersort_infiltration.png (图1: 浸润热图)\n"))
cat(sprintf("  figures/heatmap_cibersort_correlation.png (图2: 22细胞相关热图)\n"))
cat(sprintf("  figures/violin_cibersort_comparison.png (图3: 组间小提琴图)\n"))
cat(sprintf("  figures/heatmap_hub_immune_cibersort.png (图4: Hub-免疫相关热图)\n"))
cat(sprintf("  figures/pca_cibersort.png (图5: PCA聚类图)\n"))
