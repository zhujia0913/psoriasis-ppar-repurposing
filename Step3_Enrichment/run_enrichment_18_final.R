#!/usr/bin/env Rscript
# ============================================================
# 18 Human PPARgene ∩ DEG 富集分析 — 4格气泡图
# 参考图风格：横轴GeneRatio, 点大小Count, 颜色p.adjust
# 输出：go_kegg_dotplot_4panel.png
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(httr)

# 延长 KEGG 超时 (base R 方法)
options(timeout = 300)
options(download.file.method = "libcurl")

# ── Resolve script directory (portable, repo-relative) ──
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# ── 18 基因 ────────────────────────────────────────────────
genes_18 <- c("ADIPOQ", "ALDH1A3", "ANGPTL4", "APOE", "BCL2",
              "CXCL13", "CXCR4", "FABP5", "FADS2", "HBEGF",
              "HMOX1", "INSIG1", "LDLR", "MOGAT1", "NAMPT",
              "PDK4", "PLIN1", "RBP4")

# Symbol → Entrez
entrez_ids <- bitr(genes_18, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
cat("Mapped:", nrow(entrez_ids), "/", length(genes_18), "genes to Entrez\n")
if (nrow(entrez_ids) < length(genes_18)) {
  cat("Unmapped:", setdiff(genes_18, entrez_ids$SYMBOL), "\n")
}
eg <- entrez_ids$ENTREZID

# ── GO 富集 (p<0.05, q<0.05) ──────────────────────────────
go_bp <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)
go_cc <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "CC",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)
go_mf <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "MF",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05)

cat("\nGO Results:\n")
cat("  BP:", if(!is.null(go_bp)) nrow(go_bp) else 0, "terms\n")
cat("  CC:", if(!is.null(go_cc)) nrow(go_cc) else 0, "terms\n")
cat("  MF:", if(!is.null(go_mf)) nrow(go_mf) else 0, "terms\n")

# ── KEGG 富集 (p<0.05, with retry) ───────────────────────
get_kegg <- function(eg, max_retries = 3, wait_sec = 15) {
  for (i in 1:max_retries) {
    cat("  Attempt", i, "/", max_retries, "...\n")
    res <- tryCatch(
      enrichKEGG(gene = eg, organism = "hsa",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05),
      error = function(e) { cat("  Failed:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(res)) return(res)
    if (i < max_retries) {
      cat("  Retrying in", wait_sec, "s...\n")
      Sys.sleep(wait_sec)
    }
  }
  return(NULL)
}
go_kegg <- get_kegg(eg)

cat("  KEGG:", if(!is.null(go_kegg)) nrow(go_kegg) else 0, "terms\n")

# ── 保存 TSV ──────────────────────────────────────────────
outdir <- file.path(here, "data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.table(as.data.frame(go_bp),  file.path(outdir, "GO_BP_enrichment.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(go_cc),  file.path(outdir, "GO_CC_enrichment.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(go_mf),  file.path(outdir, "GO_MF_enrichment.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(go_kegg), file.path(outdir, "KEGG_enrichment.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# ── dotplot 辅助函数 ───────────────────────────────────────
make_dotplot <- function(enrich_obj, title, max_terms = 10, 
                         low_color = "#00468B", high_color = "#ED0000") {
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    return(ggplot() + annotate("text", x=0.5, y=0.5, label="No enriched terms") +
           theme_void())
  }
  n <- min(max_terms, nrow(enrich_obj))
  df <- enrich_obj@result[1:n, ]
  # 计算 GeneRatio 数值
  df$GeneRatio_num <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df <- df[order(df$p.adjust, decreasing = TRUE), ]  # reverse for top-to-bottom
  
  ggplot(df, aes(x = GeneRatio_num, y = factor(Description, levels = Description))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = low_color, high = high_color,
                         name = "p.adjust") +
    scale_size_continuous(name = "Count", range = c(3, 8)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    labs(x = "GeneRatio", y = "", title = title) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      panel.grid.major = element_line(linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

# ── 生成 4 张子图 ─────────────────────────────────────────
p_bp   <- make_dotplot(go_bp,   "(a) GO Biological Process", max_terms = 10)
p_cc   <- make_dotplot(go_cc,   "(b) GO Cellular Component", max_terms = 10)
p_mf   <- make_dotplot(go_mf,   "(c) GO Molecular Function", max_terms = 10)
p_kegg <- make_dotplot(go_kegg, "(d) KEGG Pathway",          max_terms = 10)

# ── 拼接 2×2 + 保存 ───────────────────────────────────────
library(patchwork)
combined <- (p_bp | p_cc) / (p_mf | p_kegg) +
  plot_annotation(
    title = "Enrichment Analysis of 18 Human PPARgene-DEG Overlap Genes",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(outdir, "go_kegg_dotplot_4panel.png"),
       plot = combined, width = 14, height = 10, dpi = 300, bg = "white")
cat("\n✅ 4-panel dotplot saved to:", file.path(outdir, "go_kegg_dotplot_4panel.png"), "\n")
