#!/usr/bin/env Rscript
# 重绘 Hub Genes × Immune Cells 热图
# 修正：聚类树 + Rho数值+星号 + 自适应字体颜色 + FDR标注 + 300 dpi
# 方案：pheatmap silent+main参数 + 修改text grob的gp$col颜色向量

library(pheatmap)
library(reshape2)
library(dplyr)
library(grid)
library(gtable)

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# === 1. 读取数据 ===
corr_data <- read.delim(file.path(here, "results", "hub_immune_correlation_cibersort.tsv"))

cor_wide <- dcast(corr_data, Gene ~ CellType, value.var = "Rho")
rownames(cor_wide) <- cor_wide$Gene
cor_wide <- cor_wide[, -1]
cor_matrix <- as.matrix(cor_wide)

pval_wide <- dcast(corr_data, Gene ~ CellType, value.var = "P.adj")
rownames(pval_wide) <- pval_wide$Gene
pval_wide <- pval_wide[, -1]
pval_matrix <- as.matrix(pval_wide)

# === 2. 生成显示文字 ===
display_nums <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix),
                       dimnames = dimnames(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    rho_str <- sprintf("%.2f", cor_matrix[i, j])
    star <- ifelse(pval_matrix[i, j] < 0.001, "***",
                   ifelse(pval_matrix[i, j] < 0.01, "**",
                          ifelse(pval_matrix[i, j] < 0.05, "*", "")))
    display_nums[i, j] <- paste0(rho_str, star)
  }
}

# === 3. pheatmap silent模式 + main参数 ===
ph <- pheatmap(
  cor_matrix,
  display_numbers = display_nums,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  treeheight_row = 50,
  treeheight_col = 50,
  fontsize_number = 7,
  number_color = "black",
  color = colorRampPalette(c("#2166AC", "#92B4D8", "#F7F7F7", "#D88080", "#B2182B"))(100),
  breaks = seq(-0.6, 0.6, length.out = 101),
  main = "Hub Genes x Immune Cells Spearman Correlation (CIBERSORT)\nNumbers: Rho; * FDR<0.05, ** FDR<0.01, *** FDR<0.001 (BH-FDR)",
  silent = TRUE
)

gt <- ph$gtable

# === 4. 修改text grob的字体颜色（自适应） ===
mat_idx <- which(gt$layout$name == "matrix")
mat_grob <- gt$grobs[[mat_idx]]

text_child_idx <- which(sapply(mat_grob$children, inherits, "text"))
if (length(text_child_idx) == 0) stop("No text grob found in matrix")
text_child <- mat_grob$children[[text_child_idx]]

n_row <- nrow(cor_matrix)
n_labels <- length(text_child$label)

# 基于颜色映射确定文字颜色
# pheatmap的breaks从-0.6到0.6，100色
# RdBu色板：深蓝(负值) -> 浅灰(0) -> 深红(正值)
# |Rho| >= 0.30 时背景足够深，用白字
RHO_COLOR_THRESHOLD <- 0.30

text_colors <- rep("black", n_labels)
for (k in 1:n_labels) {
  col_idx <- ((k - 1) %/% n_row) + 1
  row_from_bottom <- ((k - 1) %% n_row) + 1
  row_idx <- n_row - row_from_bottom + 1
  if (abs(cor_matrix[row_idx, col_idx]) >= RHO_COLOR_THRESHOLD) {
    text_colors[k] <- "white"
  }
}

cat(sprintf("Threshold: >= %.2f, White text: %d / %d cells\n",
            RHO_COLOR_THRESHOLD, sum(text_colors == "white"), length(text_colors)))

text_child$gp$col <- text_colors
mat_grob$children[[text_child_idx]] <- text_child
gt$grobs[[mat_idx]] <- mat_grob

# === 5. 保存 ===
outpath <- file.path(here, "figures", "heatmap_hub_immune_cibersort.png")
dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)

png(outpath, width = 14, height = 8, units = "in", res = 300)
grid.draw(gt)
dev.off()

cat("Saved:", outpath, "\n")
