#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# KEGG Enrichment — 1×2 (bar + bubble), style matched to GO 3×2
# ═══════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)
library(gtable)
library(showtext)
library(scales)
library(stringr)
showtext_auto()

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# ── Data ──────────────────────────────────────────────────
data_dir <- file.path(here, "data")
df <- read.delim(file.path(data_dir, "KEGG_enrichment.tsv"), stringsAsFactors = FALSE)

df$GR <- sapply(strsplit(df$GeneRatio, "/"),
                function(x) as.numeric(x[1]) / as.numeric(x[2]))
df <- df %>% arrange(p.adjust) %>%
  mutate(Desc_wrap = factor(str_wrap(Description, 50),
                            levels = rev(str_wrap(Description, 50))))

# ── Colors (same as GO) ───────────────────────────────────
lo_blue <- "#2B5C8F"
hi_red  <- "#D73027"

# ── Limits ────────────────────────────────────────────────
bar_xmax  <- max(df$Count) + 0.8
bub_xmax  <- max(df$GR) * 1.15
size_lims <- range(df$Count)
pval_lims <- range(df$pvalue)  # match aes(color=pvalue)
pval_breaks <- c(1e-07, 1e-06, 1e-05, 1e-04)  # ticks within pvalue range

# ── Themes (same as GO) ───────────────────────────────────
th_bar <- theme_bw(9) + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 8, color = "black"),
  axis.title.y = element_blank(),
  legend.position = "none",
  plot.margin = margin(4, 6, 2, 2)
)

th_bub <- theme_bw(9) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(linewidth = 0.3, color = "grey80"),
  axis.text = element_text(size = 8, color = "black"),
  axis.title.y = element_blank(),
  legend.position = "none",
  plot.margin = margin(4, 6, 2, 2)
)

lg_th <- theme(
  legend.position = "left",
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  legend.key.size = unit(3, "mm"),
  legend.spacing.y = unit(1, "mm"),
  legend.background = element_blank(),
  legend.margin = margin(0, 0, 0, 0)
)

# ── Plots ─────────────────────────────────────────────────
p_bar <- ggplot(df, aes(Count, Desc_wrap, fill = pvalue)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = hi_red, high = lo_blue, limits = pval_lims, trans = "log10", breaks = pval_breaks) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, bar_xmax)) +
  labs(x = "Count") + th_bar

p_bub <- ggplot(df, aes(GR, Desc_wrap, size = Count, color = pvalue)) +
  geom_point() +
  scale_color_gradient(low = hi_red, high = lo_blue, limits = pval_lims, trans = "log10", breaks = pval_breaks) +
  scale_size_continuous(range = c(2.5, 7), limits = size_lims, breaks = sort(unique(df$Count))) +
  scale_x_continuous(limits = c(0, bub_xmax),
                     labels = label_number(0.01)) +
  labs(x = "GeneRatio") + th_bub

# ── Legend: split pvalue (left-bottom) and Count (right) ──
# pvalue colorbar legend
p_leg_pval <- p_bub +
  guides(
    color = guide_colorbar(barwidth = 0.5, barheight = 3,
                           title = "pvalue", title.position = "top",
                           label.theme = element_text(size = 7)),
    size  = "none"
  ) + lg_th
leg_pval <- get_legend(p_leg_pval)

# Count size legend (keep on right)
p_leg_cnt <- p_bub +
  guides(
    color = "none",
    size  = guide_legend(title = "Count", title.position = "top")
  ) + lg_th
leg_cnt <- get_legend(p_leg_cnt)

# ── Assemble: bar + bubble + legends on right side ──────
row1 <- plot_grid(p_bar, p_bub, ncol = 2, rel_widths = c(1, 1))

final <- ggdraw() +
  draw_plot(row1, x = 0, y = 0, width = 0.82, height = 1) +
  # pvalue colorbar → right side, top (vertically aligned center x=0.91)
  draw_plot(leg_pval, x = 0.85, y = 0.55, width = 0.12, height = 0.40) +
  # Count size legend → right side, bottom (same center x as pvalue)
  draw_plot(leg_cnt, x = 0.85, y = 0.10, width = 0.12, height = 0.40) +

  # Panel letters
  draw_label("a", x = 0.02, y = 0.96, size = 15, fontface = "bold", hjust = 0) +
  draw_label("b", x = 0.44, y = 0.96, size = 15, fontface = "bold", hjust = 0)

# ── Save ──────────────────────────────────────────────────
out <- file.path(here, "figures", "Figure_KEGG_Enrichment")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(out, ".png"), final, width = 14, height = 4, dpi = 300, bg = "white")
pdf(paste0(out, ".pdf"), width = 14, height = 4, bg = "white"); print(final); dev.off()

cat("\n✅ Done\n")
system(sprintf("sips -g pixelHeight -g pixelWidth '%s.png'", out))
cat(sprintf("  PNG: %s.png\n  PDF: %s.pdf\n", out, out))
