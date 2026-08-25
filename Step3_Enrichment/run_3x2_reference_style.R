#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════
# GO Enrichment 3×2 v6 — row-wise cowplot + single legend stack
# Font: Arial (all text elements)
# ═══════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)
library(gtable)
library(showtext)
showtext_auto()

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# ── Register Arial font ───────────────────────────────────
font_add("Arial",
         regular    = "/System/Library/Fonts/Supplemental/Arial.ttf",
         bold       = "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         italic     = "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")

# ── Data ──────────────────────────────────────────────────
data_dir <- file.path(here, "data")
read_top10 <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE)
  df$GR <- sapply(strsplit(df$GeneRatio, "/"),
                   function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df <- df %>% arrange(p.adjust) %>% head(10) %>%
    mutate(Desc_wrap = factor(stringr::str_wrap(Description, 50),
                              levels = rev(stringr::str_wrap(Description, 50))))
  df
}

bp <- read_top10(file.path(data_dir, "GO_BP_enrichment.tsv"))
cc <- read_top10(file.path(data_dir, "GO_CC_enrichment.tsv"))
mf <- read_top10(file.path(data_dir, "GO_MF_enrichment.tsv"))

lo_blue  <- "#2B5C8F"   # high pvalue
hi_red   <- "#D73027"   # low pvalue

# ── Global limits ─────────────────────────────────────────
bar_xmax  <- max(c(bp$Count, cc$Count, mf$Count)) + 0.8
bub_xmax  <- max(c(bp$GR, cc$GR, mf$GR)) * 1.08
size_lims <- range(c(bp$Count, cc$Count, mf$Count))

# ── Themes ────────────────────────────────────────────────
th_bar <- theme_bw(9, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size=8, color="black", family="Arial"),
  axis.title.y = element_blank(),
  legend.position = "none",
  plot.margin = margin(4, 6, 2, 2)
)

th_bub <- theme_bw(9, base_family = "Arial") + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(linewidth=0.3, color="grey80"),
  axis.text = element_text(size=8, color="black", family="Arial"),
  axis.title.y = element_blank(),
  legend.position = "none",
  plot.margin = margin(4, 6, 2, 2)
)

lg_th <- theme(
  legend.position = "right",
  legend.title = element_text(size=8, family="Arial"),
  legend.text = element_text(size=7, family="Arial"),
  legend.key.size = unit(3, "mm"),
  legend.spacing.y = unit(1, "mm"),
  legend.background = element_blank(),
  legend.margin = margin(0,0,0,0)
)

# ── Plot builders ─────────────────────────────────────────
bar_plot <- function(df) {
  ggplot(df, aes(Count, Desc_wrap, fill=p.adjust)) +
    geom_col(width=0.7) +
    scale_fill_gradient(low=hi_red, high=lo_blue) +
    scale_x_continuous(expand=expansion(mult=c(0, 0.06)), limits=c(0, bar_xmax)) +
    labs(x="Count") + th_bar
}

bub_plot <- function(df) {
  ggplot(df, aes(GR, Desc_wrap, size=Count, color=p.adjust)) +
    geom_point() +
    scale_color_gradient(low=hi_red, high=lo_blue) +
    scale_size_continuous(range=c(2.5, 7), limits=size_lims) +
    scale_x_continuous(limits=c(0, bub_xmax),
                       labels=scales::label_number(0.01)) +
    labs(x="GeneRatio") + th_bub
}

# 6 subplots (NO legends)
bar_bp  <- bar_plot(bp);  bub_bp  <- bub_plot(bp)
bar_cc  <- bar_plot(cc);  bub_cc  <- bub_plot(cc)
bar_mf  <- bar_plot(mf);  bub_mf  <- bub_plot(mf)

# ── Legend: combined Count + p.adjust, vertically stacked & aligned ──
pval_all_lims <- range(c(bp$p.adjust, cc$p.adjust, mf$p.adjust))

# Single plot holding BOTH guides -> ggplot2 stacks them vertically,
# guarantees horizontal alignment (same legend box, centered)
bub_leg_combined <- bub_cc +
  scale_color_gradient(low=hi_red, high=lo_blue, limits=pval_all_lims) +
  scale_size_continuous(range=c(2.5, 7), limits=size_lims) +
  guides(
    size  = guide_legend(title="Count", title.position="top",
                         order=1),
    color = guide_colorbar(barwidth=0.4, barheight=6,
                           title="p.adjust", title.position="top",
                           label.theme=element_text(size=7, family="Arial"),
                           direction="vertical",
                           title.theme=element_text(size=8, hjust=0.5, family="Arial"),
                           order=2)
  ) +
  theme_void() +
  theme(legend.position="right",
        legend.title=element_text(size=8, family="Arial"),
        legend.text=element_text(size=7, family="Arial"),
        legend.key.size=unit(3,"mm"),
        legend.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0,"mm"),
        legend.box="vertical",
        legend.box.margin=margin(0,0,0,0))
leg_combined <- get_legend(bub_leg_combined)

# ── Align left/right columns independently for uniform widths ──
bars <- align_plots(bar_bp, bar_cc, bar_mf, align="v", axis="lr")
bubs <- align_plots(bub_bp, bub_cc, bub_mf, align="v", axis="lr")

# ── Assemble rows (no per-row align; widths already matched) ──
row1 <- plot_grid(bars[[1]], bubs[[1]], ncol=2, rel_widths=c(1,1))
row2 <- plot_grid(bars[[2]], bubs[[2]], ncol=2, rel_widths=c(1,1))
row3 <- plot_grid(bars[[3]], bubs[[3]], ncol=2, rel_widths=c(1,1))

grid_all <- plot_grid(row1, row2, row3, ncol=1, align="v", axis="lr")

# ── Final: combined legend stack on the right side ─
# Use grid graphics for precise legend placement
final <- ggdraw() +
  # Main 6-panel grid (uses left ~82% of width, leave room for legends)
  draw_plot(grid_all, x=0, y=0, width=0.82, height=1) +
  # Combined legend (Count on top, p.adjust below) — vertically aligned
  draw_plot(leg_combined, x=0.86, y=0.18, width=0.13, height=0.64) +
  # Ontology labels -- outside bubble panels, left of legend
  draw_label("BP", x=0.855, y=0.833, size=12, fontface="bold", hjust=0.5, fontfamily="Arial") +
  draw_label("CC", x=0.855, y=0.500, size=12, fontface="bold", hjust=0.5, fontfamily="Arial") +
  draw_label("MF", x=0.855, y=0.167, size=12, fontface="bold", hjust=0.5, fontfamily="Arial") +
  # Panel letters
  draw_label("a", x=0.02, y=0.985, size=15, fontface="bold", hjust=0, fontfamily="Arial") +
  draw_label("b", x=0.425, y=0.985, size=15, fontface="bold", hjust=0, fontfamily="Arial")

# ── Save ─────────────────
out <- file.path(here, "figures", "Figure_GO_Enrichment_3x2")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

# Save PNG @ 3000 DPI (14in x 3000 = 42000 px wide)
# (R on macOS quartz does not support LZW for TIFF, so PNG is produced here
#  and then converted to LZW-compressed TIFF via fix_dpi_tiff.py)
png(paste0(out, ".png"), width=14, height=12, units="in", res=3000, bg="white")
print(final)
dev.off()

# Save PDF (vector, resolution-independent; showtext renders text as paths)
pdf(paste0(out, ".pdf"), width=14, height=12, bg="white")
print(final)
dev.off()

# Post-process: embed 3000 DPI metadata into PNG + write LZW TIFF
system(sprintf("python3 '%s/fix_dpi_tiff.py'", here))

cat("\n✅ Done\n")
cat(sprintf("  PNG:  %s.png  (3000 DPI)\n", out))
cat(sprintf("  TIFF: %s.tiff (3000 DPI, LZW)\n", out))
cat(sprintf("  PDF:  %s.pdf  (vector)\n", out))
