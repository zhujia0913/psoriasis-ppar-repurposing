# =============================================================================
# CIBERSORT 组间差异图（文献风格：散点+jitter+中位数+IQR）
# 对齐文献 Fig 6C 风格
# Font: Arial | Output: 300 DPI (BMC 单文件上限 10 MB，3000 DPI 会超限)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(grid)
  library(showtext)
})
showtext_auto()

# --- Register Arial font (all 4 variants) ---
font_add("Arial",
         regular    = "/System/Library/Fonts/Supplemental/Arial.ttf",
         bold       = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         italic     = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")

# --- 1. 路径 (BMC投递版: 相对本脚本运行位置，可在任意机器复跑) ---
resolve_base <- function() {
  for (b in c(".", "..")) {
    if (file.exists(file.path(b, "results", "cibersort_results_with_group.tsv"))) {
      return(normalizePath(b))
    }
  }
  stop("找不到 results/cibersort_results_with_group.tsv，请在包根目录或 scripts/ 下运行本脚本")
}
base_dir <- resolve_base()
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
scripts_dir <- file.path(base_dir, "scripts")

# --- 2. 读取数据（保留原始列名空格） ---
results_with_group <- read.delim(file.path(results_dir, "cibersort_results_with_group.tsv"),
                                  check.names = FALSE)
diff_df <- read.delim(file.path(results_dir, "cibersort_comparison_NN_vs_PP.tsv"),
                      check.names = FALSE)

# --- 3. 细胞类型（按显著性排序） ---
cell_cols <- diff_df$CellType

# --- 4. 计算每个细胞的 y.position（顶部括号位置） ---
y_pos_df <- results_with_group %>%
  select(all_of(cell_cols)) %>%
  summarise(across(everything(), ~ max(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "CellType", values_to = "max_val") %>%
  mutate(y.position = max_val + 0.05)

# --- 5. 构造 sig_df ---
sig_df <- diff_df %>%
  mutate(
    label = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01  ~ "**",
      P.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    group1 = "NN",
    group2 = "PP",
    CellType = factor(CellType, levels = cell_cols)
  ) %>%
  left_join(y_pos_df %>% select(CellType, y.position), by = "CellType")

# --- 6. 长格式数据 ---
results_long <- results_with_group %>%
  select(Sample, Group, all_of(cell_cols)) %>%
  pivot_longer(cols = all_of(cell_cols),
               names_to = "CellType",
               values_to = "Proportion")

results_long$CellType <- factor(results_long$CellType, levels = cell_cols)
results_long$Group    <- factor(results_long$Group,    levels = c("NN", "PP"))

# --- 7. 绘图 ---
p <- ggplot(results_long, aes(x = CellType, y = Proportion, color = Group)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6),
              size = 0.8, alpha = 0.65) +
  stat_summary(fun = median,
               fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75),
               geom = "errorbar", width = 0.25, linewidth = 0.6,
               position = position_dodge(width = 0.6), color = "black") +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6,
               position = position_dodge(width = 0.6), color = "black") +
  geom_text(data = sig_df,
            aes(x = CellType, y = y.position, label = label),
            inherit.aes = FALSE, size = 3.5, fontface = "bold",
            color = "black", vjust = 0, family = "Arial") +
  scale_color_manual(values = c("NN" = "#4DBBD5", "PP" = "#E64B35")) +
  scale_y_continuous(limits = c(0, max(y_pos_df$y.position) + 0.12),
                     breaks = pretty(c(0, max(y_pos_df$y.position) + 0.12), 7)) +
  labs(x = NULL, y = "Fraction", color = NULL) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 9, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 11, color = "black", family = "Arial"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    plot.margin = margin(10, 15, 10, 10)
  )

# --- 右上角图例：NN/PP（左）+ FDR 阈值（右）---  "Wilcoxon + BH-FDR" 已删除 ---
legend_grob <- grobTree(
  rectGrob(x = 0, y = 0, width = 1, height = 1,
           just = c("left", "bottom"),
           gp = gpar(col = "grey70", fill = "white", lwd = 0.3)),
  # 左列：NN/PP
  pointsGrob(x = 0.08, y = 0.85, pch = 21, size = unit(0.25, "cm"),
             gp = gpar(fill = "#4DBBD5", col = "#4DBBD5")),
  textGrob("NN (n=91)", x = 0.18, y = 0.87, just = c("left", "top"),
           gp = gpar(fontsize = 7.5, fontfamily = "Arial")),
  pointsGrob(x = 0.08, y = 0.62, pch = 21, size = unit(0.25, "cm"),
             gp = gpar(fill = "#E64B35", col = "#E64B35")),
  textGrob("PP (n=118)", x = 0.18, y = 0.64, just = c("left", "top"),
           gp = gpar(fontsize = 7.5, fontfamily = "Arial")),
  # 右列：FDR 阈值（底部留空以替代原先的"Wilcoxon + BH-FDR"行）
  textGrob("*  FDR < 0.05",    x = 0.50, y = 0.87, just = c("left", "top"),
           gp = gpar(fontsize = 7.5, fontfamily = "Arial")),
  textGrob("**  FDR < 0.01",   x = 0.50, y = 0.64, just = c("left", "top"),
           gp = gpar(fontsize = 7.5, fontfamily = "Arial")),
  textGrob("***  FDR < 0.001", x = 0.50, y = 0.41, just = c("left", "top"),
           gp = gpar(fontsize = 7.5, fontfamily = "Arial"))
)

p <- p + annotation_custom(legend_grob,
                           xmin = 19.0, xmax = 22,
                           ymin = 0.50, ymax = 0.72)

# --- 8. 保存 ---
# 物理尺寸 16 x 6 inches (与 PDF 一致)
# BMC 要求约 300 dpi 即可；3000 dpi 会产生 >10 MB 的单文件，被投稿系统拒绝

# PNG @ 300 DPI -> 4800 x 1800 pixels (远低于 10 MB，且满足 300 dpi)
png(file.path(figures_dir, "violin_cibersort_19cells_literature.png"),
    width = 16, height = 6, units = "in", res = 300, bg = "white")
print(p)
dev.off()

# TIFF @ 300 DPI with LZW compression (R 原生支持)
tiff(file.path(figures_dir, "violin_cibersort_19cells_literature.tiff"),
     width = 16, height = 6, units = "in", res = 300, compression = "lzw", bg = "white")
print(p)
dev.off()

# PDF (vector, resolution-independent)
pdf(file.path(figures_dir, "violin_cibersort_19cells_literature.pdf"),
    width = 16, height = 6, bg = "white")
print(p)
dev.off()

cat("\n✅ Done\n")
cat(sprintf("  PNG:  %s/violin_cibersort_19cells_literature.png  (300 DPI)\n", figures_dir))
cat(sprintf("  TIFF: %s/violin_cibersort_19cells_literature.tiff (300 DPI, LZW)\n", figures_dir))
cat(sprintf("  PDF:  %s/violin_cibersort_19cells_literature.pdf  (vector)\n", figures_dir))
cat(sprintf("  细胞数: %d\n", length(cell_cols)))
cat("  显著标注:", paste(sig_df$label, collapse = " "), "\n")
