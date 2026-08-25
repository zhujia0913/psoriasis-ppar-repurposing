#!/usr/bin/env Rscript
# =============================================================================
# 银屑病 DEG 火山图 — GSE13355/GSE14905/GSE78097 三数据集联合分析
# 样式对齐用户提供的参考图
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

set.seed(42)

# ---- Resolve script directory (portable, repo-relative) ----
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# ---- 参数 ----
DEG_FILE  <- file.path(here, "deg_results_all_genes.tsv")
OUT_DIR   <- file.path(here, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 读取 ----
deg <- read.table(DEG_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Total genes: %d\n", nrow(deg)))

# ---- 阈值定义 ----
logFC_thr <- 1
fdr_thr   <- 0.05

deg <- deg %>%
  mutate(
    neg_log10_fdr = -log10(adj.P.Val),
    category = case_when(
      adj.P.Val < fdr_thr & logFC >= logFC_thr  ~ "Up",
      adj.P.Val < fdr_thr & logFC <= -logFC_thr ~ "Down",
      TRUE ~ "NS"
    )
  )

cat(sprintf("Up: %d, Down: %d, NS: %d\n",
            sum(deg$category == "Up"),
            sum(deg$category == "Down"),
            sum(deg$category == "NS")))

# ---- 标注基因 ----
# 取 |logFC| 最大的上调和下调各 top N
N_label_each <- 10
up_deg   <- deg %>% filter(category == "Up") %>% arrange(desc(logFC))
down_deg <- deg %>% filter(category == "Down") %>% arrange(logFC)

label_genes <- bind_rows(
  head(up_deg, N_label_each),
  head(down_deg, N_label_each)
)

cat(sprintf("Labeled genes: %d\n", nrow(label_genes)))
cat("Up labels:", paste(head(up_deg$Gene, N_label_each), collapse = ", "), "\n")
cat("Down labels:", paste(head(down_deg$Gene, N_label_each), collapse = ", "), "\n")

# ---- 配色 ----
# 参考图风格：灰/红/蓝
colors <- c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "#B0B0B0")

# ---- Plotting ----
# 限制 y 轴裁剪极端值以提高可读性，但不影响 p 值
y_cap <- max(deg$neg_log10_fdr[is.finite(deg$neg_log10_fdr)]) * 1.02
x_lim <- c(-max(abs(deg$logFC)) * 1.08, max(abs(deg$logFC)) * 1.08)

cat(sprintf("Y cap: %.1f, X lim: ±%.1f\n", y_cap, max(x_lim)))

# 打乱绘制顺序避免点重叠遮盖
deg_plot <- deg[sample(nrow(deg)), ]

p <- ggplot(deg_plot, aes(x = logFC, y = neg_log10_fdr)) +
  # NS 点先画
  geom_point(data = subset(deg_plot, category == "NS"),
             aes(color = category), size = 1.2, alpha = 0.50) +
  # 下调点
  geom_point(data = subset(deg_plot, category == "Down"),
             aes(color = category), size = 1.5, alpha = 0.75) +
  # 上调点
  geom_point(data = subset(deg_plot, category == "Up"),
             aes(color = category), size = 1.5, alpha = 0.75) +
  
  # 阈值线
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed",
             color = "#333333", linewidth = 0.4) +
  geom_vline(xintercept = c(-logFC_thr, logFC_thr), linetype = "dashed",
             color = "#333333", linewidth = 0.4) +
  
  # 基因标签（已删除）
  # geom_text_repel(
  #   data = label_genes,
  #   aes(label = Gene, color = category),
  #   size = 3.2,
  #   max.overlaps = 25,
  #   box.padding = 0.35,
  #   point.padding = 0.25,
  #   segment.size = 0.25,
  #   segment.alpha = 0.6,
  #   show.legend = FALSE,
  #   fontface = "italic"
  # ) +
  
  # 颜色
  scale_color_manual(values = colors, name = "",
                     breaks = c("Up", "Down", "NS")) +
  
  # 坐标轴
  scale_x_continuous(
    limits = x_lim,
    breaks = seq(-8, 9, by = 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, y_cap),
    expand = c(0, 0)
  ) +
  
  # 标签
  labs(
    x = expression(log[2] * "(Fold Change)"),
    y = expression(-log[10] * "(FDR)"),
    title = "Psoriasis DEGs (GSE13355 + GSE14905 + GSE78097)"
  ) +
  
  # 主题
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain",
                              margin = margin(b = 12)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "grey70",
                                     linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10),
    legend.margin = margin(4, 6, 4, 6),
    plot.margin = margin(15, 15, 10, 10),
    aspect.ratio = 0.85
  ) +

  # 左上角颜色图例
  guides(color = guide_legend(
    title = NULL,
    override.aes = list(size = 3, alpha = 1)
  ))

# ---- 保存 ----
ggsave(file.path(OUT_DIR, "Volcano_DEG_three_datasets.png"),
       plot = p, width = 9, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "Volcano_DEG_three_datasets.pdf"),
       plot = p, width = 9, height = 8.5, bg = "white")

cat(sprintf("\nDone — saved to %s\n", OUT_DIR))
