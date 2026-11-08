# =============================================================================
# Step 9 — Section 11 only: Hub gene x Immune cell Spearman correlation (rerun)
# 使用最终 CIBERSORT 比例 (cibersort_proportions_all.tsv, 2026-11-02 重新生成)
# + 最终 Pre-ComBat 表达矩阵，重算 4 hub genes x 22 免疫细胞 相关并做 BH-FDR。
# 修正陈旧的 hub_immune_correlation_cibersort.tsv (原 Jun 29 版本不一致)。
# =============================================================================

# --- 路径设置 (自动适配运行目录, 可在任意机器直接复跑) ---
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

library(pheatmap)

# --- 读取最终表达矩阵 (基因行, 样本列) ---
expr_pre <- read.csv(file.path(base_dir, "raw_data", "expression_pre_combat.csv"),
                     row.names = 1, check.names = FALSE)
cat(sprintf("Pre-ComBat表达矩阵: %d基因 x %d样本\n", nrow(expr_pre), ncol(expr_pre)))

# --- 读取最终 CIBERSORT 比例 (209 样本 x 22 细胞) ---
results_all <- as.matrix(read.delim(file.path(results_dir, "cibersort_proportions_all.tsv"),
                                    row.names = 1, check.names = FALSE))
cell_cols <- colnames(results_all)
cat(sprintf("CIBERSORT比例: %d样本 x %d细胞\n", nrow(results_all), ncol(results_all)))

# --- 对齐样本顺序 (比例行序 = CIBERSORT 输出序 = common_samples 序) ---
common_samples <- intersect(colnames(expr_pre), rownames(results_all))
expr_pre <- expr_pre[, common_samples, drop = FALSE]
results_all <- results_all[common_samples, , drop = FALSE]
cat(sprintf("对齐后共有样本: %d\n", length(common_samples)))

# --- Hub 基因 ---
hub_genes <- c("ADIPOQ", "APOE", "PLIN1", "BCL2")
hub_expr <- expr_pre[hub_genes, , drop = FALSE]

# --- Spearman 相关 + p ---
# 注意: R 的 cor.test(method="spearman") 在 n=209 且相关较强时会用精确法计算 p 值并
# 数值下溢为 0 (例如 PLIN1-Mast cells resting 真实 p≈5.5e-13 被记为 0)。
# 论文中 p=0 / FDR=0 不规范, 故改用与大样本 t 近似 (scipy spearmanr 同法) 计算 p:
#   t = rho * sqrt((n-2)/(1-rho^2));  p = 2 * pt(-|t|, n-2)
cor_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols),
                  dimnames = list(hub_genes, cell_cols))
p_hub <- matrix(NA, nrow = length(hub_genes), ncol = length(cell_cols),
                dimnames = list(hub_genes, cell_cols))

n_samp <- nrow(results_all)
for (gene in hub_genes) {
  for (ct in cell_cols) {
    gene_vec <- as.numeric(hub_expr[gene, ])
    rho <- cor(gene_vec, results_all[, ct], method = "spearman")
    tstat <- rho * sqrt((n_samp - 2) / (1 - rho^2))
    pval <- 2 * pt(-abs(tstat), n_samp - 2)
    cor_hub[gene, ct] <- rho
    p_hub[gene, ct] <- pval
  }
}

# --- BH-FDR 校正 ---
p_hub_adj <- matrix(p.adjust(as.vector(p_hub), method = "BH"),
                    nrow = length(hub_genes), ncol = length(cell_cols),
                    dimnames = list(hub_genes, cell_cols))

# --- 保存 ---
# 注意: cor_hub/p_hub/p_hub_adj 为 行=hub基因, 列=免疫细胞 的矩阵。
# R 的 as.vector() 按列展开, 而 Gene/CellType 按行构造, 直接 as.vector() 会错位。
# 必须用 t() 转置后按行展开, 才能保证 Rho/P 与标签一一对应。
hub_cor_df <- data.frame(
  Gene = rep(hub_genes, each = length(cell_cols)),
  CellType = rep(cell_cols, length(hub_genes)),
  Rho = as.vector(t(cor_hub)),
  P.value = as.vector(t(p_hub)),
  P.adj = as.vector(t(p_hub_adj))
)
hub_cor_df <- hub_cor_df[order(hub_cor_df$P.adj), ]
out_tsv <- file.path(results_dir, "hub_immune_correlation_cibersort.tsv")
write.table(hub_cor_df, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("修正版 hub-免疫相关保存: %s\n", out_tsv))

# --- Fig 13: Hub-免疫相关热图 ---
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
  main = "Hub Genes x Immune Cells Spearman Correlation (CIBERSORT)"
)
dev.off()
cat("图4(Fig13)保存: heatmap_hub_immune_cibersort.png\n")

pdf(file.path(figures_dir, "heatmap_hub_immune_cibersort.pdf"),
    width = 12, height = 6)
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
  main = "Hub Genes x Immune Cells Spearman Correlation (CIBERSORT)"
)
dev.off()
cat("图4(Fig13)保存: heatmap_hub_immune_cibersort.pdf\n")

# --- 汇总 ---
cat(sprintf("\nHub基因显著相关 (P.adj < 0.05): %d / %d\n",
            sum(p_hub_adj < 0.05), length(p_hub_adj)))
cat("BCL2-Tregs: Rho=", round(cor_hub["BCL2","T cells regulatory (Tregs)"],4),
    " P.adj=", format.pval(p_hub_adj["BCL2","T cells regulatory (Tregs)"],4), "\n")
