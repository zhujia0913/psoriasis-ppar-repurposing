# 仅重画 Step9 图件 (Fig3 小提琴 / Fig4 Hub-免疫热图)，不重跑 CIBERSORT。
# 从最终数据读取，使用与 run_cibersort_analysis_rerun.R 一致且已修复的算法。
# 仅重画 Step9 图件 (Fig3 小提琴 / Fig4 Hub-免疫热图)，不重跑 CIBERSORT。
# 从最终数据读取，使用与 run_cibersort_analysis_rerun.R 一致且已修复的算法。
# 路径设置 (BMC投递版: 自动适配运行目录)
resolve_base <- function() {
  for (b in c(".", "..")) {
    if (file.exists(file.path(b, "raw_data", "expression_pre_combat.csv"))) return(normalizePath(b))
  }
  stop("找不到 raw_data/expression_pre_combat.csv：请在 免疫浸润/ 根目录或其 scripts/ 子目录下运行本脚本。")
}
base_dir <- resolve_base()
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")

expr_pre <- read.csv(file.path(base_dir, "raw_data", "expression_pre_combat.csv"),
                     row.names = 1, check.names = FALSE)
sample_info <- read.delim(file.path(data_dir, "sample_info.tsv"))
names(sample_info)[1:2] <- c("Sample", "Group")
prop <- read.delim(file.path(results_dir, "cibersort_proportions_all.tsv"), row.names = 1, check.names = FALSE)
cell_cols <- colnames(prop)

# ---------- Fig3: 组间小提琴图 (用 ggplot2 geom_violin, 规避 ggpubr 新版兼容报错) ----------
library(ggplot2); library(ggpubr)
results_long <- reshape2::melt(as.matrix(prop))
colnames(results_long) <- c("Sample", "CellType", "Proportion")
results_long <- merge(results_long, sample_info[, c("Sample", "Group")], by = "Sample")
results_long$CellType <- factor(results_long$CellType, levels = cell_cols)
results_long$Group <- factor(results_long$Group, levels = c("NN", "PP"))

png(file.path(figures_dir, "violin_cibersort_comparison.png"), width = 2400, height = 1200, res = 150)
p <- ggplot(results_long, aes(x = Group, y = Proportion, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 4) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  labs(title = "Immune Cell Infiltration: NN vs PP (CIBERSORT)", x = NULL, y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 7), legend.position = "bottom")
print(p); dev.off()
cat("图3 (violin) 重画完成\n")

# ---------- Fig4: Hub-免疫 Spearman 热图 (修复版算法) ----------
hub_genes <- c("ADIPOQ", "APOE", "PLIN1", "BCL2")
hub_expr <- expr_pre[hub_genes, , drop = FALSE]
cor_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols), dimnames = list(hub_genes, cell_cols))
p_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols), dimnames = list(hub_genes, cell_cols))
n <- nrow(prop)
for (gene in hub_genes) {
  for (ct in cell_cols) {
    gene_vec <- as.numeric(hub_expr[gene, ])
    rho <- cor(gene_vec, prop[, ct], method = "spearman")
    tstat <- rho * sqrt((n - 2) / (1 - rho^2))
    pval <- 2 * pt(-abs(tstat), n - 2)
    cor_hub[gene, ct] <- rho
    p_hub[gene, ct] <- pval
  }
}
p_hub_adj <- matrix(p.adjust(as.vector(p_hub), method = "BH"),
                    nrow = length(hub_genes), ncol = length(cell_cols), dimnames = list(hub_genes, cell_cols))
sig_stars <- ifelse(p_hub_adj < 0.001, "***", ifelse(p_hub_adj < 0.01, "**", ifelse(p_hub_adj < 0.05, "*", "")))

# 图注约定：hub genes 为列(4)、immune cell types 为行(22) → 对矩阵转置
cor_plot   <- t(cor_hub)
sig_plot   <- t(sig_stars)
display_txt <- matrix(paste0(sprintf("%.2f", cor_plot), sig_plot),
                      nrow = nrow(cor_plot), ncol = ncol(cor_plot),
                      dimnames = dimnames(cor_plot))

png(file.path(figures_dir, "heatmap_hub_immune_cibersort.png"), width = 2000, height = 1300, res = 300)
pheatmap::pheatmap(cor_plot, display_numbers = display_txt, cluster_cols = TRUE, cluster_rows = TRUE,
                   fontsize_row = 12, fontsize_col = 11, fontsize_number = 9, number_color = "black",
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   breaks = seq(-0.6, 0.6, length.out = 101),
                   main = "Hub genes x Immune cells (Spearman rho)")
dev.off()

# 投稿用 TIFF (300 DPI, 单文件 < 10 MB；macOS quartz 设备不支持压缩，未压缩约 7 MB 仍合规)
tiff(file.path(figures_dir, "heatmap_hub_immune_cibersort.tiff"),
     width = 2000, height = 1300, res = 300, bg = "white")
pheatmap::pheatmap(cor_plot, display_numbers = display_txt, cluster_cols = TRUE, cluster_rows = TRUE,
                   fontsize_row = 12, fontsize_col = 11, fontsize_number = 9, number_color = "black",
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   breaks = seq(-0.6, 0.6, length.out = 101),
                   main = "Hub genes x Immune cells (Spearman rho)")
dev.off()
# 矢量 PDF (BMC 首选; 与 PNG/TIFF 同源，含 Rho 数值 + 显著性星号)
pdf(file.path(figures_dir, "heatmap_hub_immune_cibersort.pdf"), width = 10, height = 6.5)
pheatmap::pheatmap(cor_plot, display_numbers = display_txt, cluster_cols = TRUE, cluster_rows = TRUE,
                   fontsize_row = 12, fontsize_col = 11, fontsize_number = 9, number_color = "black",
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   breaks = seq(-0.6, 0.6, length.out = 101),
                   main = "Hub genes x Immune cells (Spearman rho)")
dev.off()
cat("图4 (heatmap_hub_immune) 重画完成\n")
