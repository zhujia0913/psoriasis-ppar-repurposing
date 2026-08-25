# =============================================================================
# 银屑病微阵列跨数据集整合 DEG 分析 — 标准化管线 (v2 AnnoProbe)
# 
# v2 变更: 注释源从 GEO GPL570 官方表 → AnnoProbe::idmap(type="bioc")
#   原因: GEO 旧版注释缺失 ZNF683/CAMP 等基因
#   AnnoProbe 基于 Bioconductor hgu133plus2.db 最新版 (41,937 probe-gene pairs)
# 
# 管线参照标准:
#   归一化:      RMA (Irizarry 2003, Biostatistics)
#   批次校正:    ComBat with mod (Johnson, Li & Rabinovic 2007, Biostatistics)
#   探针→基因:   AnnoProbe idmap (Zeng JM et al., v0.1.0, Bioconductor hgu133plus2.db)
#   多探针→同基因: MaxMean (Miller 2011, BMC Bioinformatics)
#   差异检验:    limma eBayes(trend=TRUE) (Ritchie 2015, NAR)
#   多重校正:    Benjamini-Hochberg FDR (B-H 1995)
#   阈值:       |logFC| > 1, FDR < 0.05
#
# 数据集: GSE13355, GSE14905, GSE78097
# =============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(sva)        # ComBat
  library(limma)      # eBayes
  library(AnnoProbe)  # v2: probe re-annotation
})

# ---- Resolve script directory (portable, repo-relative) ----
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  here <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  here <- getwd()
}

# ---- 参数 ----
# RAW_DIR: 原始 GEO CEL/RDS 文件所在目录 (GSE13355/GSE14905/GSE78097)。
# 这些 CEL 文件属于公开 GEO 下载数据，不随本仓库分发。
# 用户须将 CEL 文件放入下方默认目录，或自行设置 RAW_DIR 指向本地目录。
RAW_DIR <- file.path(here, "data", "raw")
OUT_DIR <- here

dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATASETS <- c("GSE13355", "GSE14905", "GSE78097")

# =============================================================================
# Step 1: 加载 & 筛选 NN/PP 样本
# =============================================================================
cat("=== Step 1: 加载数据 & 筛选 NN/PP ===\n")

expr_list <- list()
group_list <- list()
batch_list <- list()

for (ds in DATASETS) {
  rds_path <- file.path(DATA_DIR, paste0(ds, "_rma_only.RDS"))
  rds      <- readRDS(rds_path)
  eSet     <- rds$expr
  group    <- rds$group
  
  # 仅保留 NN 和 PP
  keep     <- group %in% c("NN", "PP")
  expr_sub <- exprs(eSet)[, keep]
  group_sub <- factor(group[keep], levels = c("NN", "PP"))
  
  expr_list[[ds]]  <- expr_sub
  group_list[[ds]] <- group_sub
  batch_list[[ds]] <- rep(ds, ncol(expr_sub))
  
  cat(sprintf("  %s: %d probes × %d samples (NN=%d, PP=%d)\n",
              ds, nrow(expr_sub), ncol(expr_sub),
              sum(group_sub == "NN"), sum(group_sub == "PP")))
}

# =============================================================================
# Step 2: 探针→基因映射 (AnnoProbe — Bioconductor hgu133plus2.db)
# =============================================================================
cat("\n=== Step 2: 探针→基因映射 (AnnoProbe) ===\n")

# v2: 使用 AnnoProbe::idmap 获取最新注释
# type="bioc" 基于 Bioconductor hgu133plus2.db 包，覆盖 41,937 probe-gene pairs
# 相比 GEO 官方表（~14,900 有效对），额外回收 ZNF683、CAMP 等在新版注释中新增的基因
annot_raw <- idmap("GPL570", type = "bioc")

# AnnoProbe 返回 clean 1-to-1 映射（无 /// 多基因探针，无需额外过滤）
probe2gene <- data.frame(
  probe = as.character(annot_raw$probe_id),
  gene  = as.character(annot_raw$symbol),
  stringsAsFactors = FALSE
)

cat(sprintf("  AnnoProbe 注释: %d 探针-基因对 (唯一基因: %d)\n",
            nrow(probe2gene), length(unique(probe2gene$gene))))

# 关键基因回收检查
cat("  关键基因回收:\n")
for (g in c("RUNX3", "ZNF683", "CAMP", "S100A7", "S100A8", "S100A9", "LCN2")) {
  hits <- probe2gene$probe[probe2gene$gene == g]
  status <- if (length(hits) > 0) paste(length(hits), "probes:", paste(hits, collapse=",")) else "MISSING"
  cat(sprintf("    %s: %s\n", g, status))
}

# =============================================================================
# Step 3: 对每个数据集: 探针替换→MaxMean 压缩
# =============================================================================
cat("\n=== Step 3: MaxMean 多探针压缩 ===\n")

collapse_probes <- function(expr_mat, p2g) {
  # 匹配探针
  matched <- p2g[p2g$probe %in% rownames(expr_mat), ]
  cat(sprintf("  匹配探针: %d / %d\n", nrow(matched), nrow(expr_mat)))
  
  expr_matched <- expr_mat[matched$probe, , drop = FALSE]
  gene_list    <- matched$gene
  
  # MaxMean: 对每个基因，取均值最高的探针
  ug <- unique(gene_list)
  gene_expr <- matrix(NA, nrow = length(ug), ncol = ncol(expr_matched))
  rownames(gene_expr) <- ug
  colnames(gene_expr) <- colnames(expr_matched)
  
  for (i in seq_along(ug)) {
    g <- ug[i]
    probes <- which(gene_list == g)
    if (length(probes) == 1) {
      gene_expr[i, ] <- expr_matched[probes, ]
    } else {
      sub_expr <- expr_matched[probes, , drop = FALSE]
      row_means <- rowMeans(sub_expr)
      best <- which.max(row_means)
      gene_expr[i, ] <- sub_expr[best, ]
    }
  }
  cat(sprintf("  → %d 基因\n", nrow(gene_expr)))
  return(gene_expr)
}

gene_expr_list <- list()
for (ds in DATASETS) {
  cat(sprintf("  %s:\n", ds))
  gene_expr_list[[ds]] <- collapse_probes(expr_list[[ds]], probe2gene)
}

# =============================================================================
# Step 4: 取基因交集 + 合并表达矩阵
# =============================================================================
cat("\n=== Step 4: 基因交集 → 合并表达矩阵 ===\n")

common_genes <- Reduce(intersect, lapply(gene_expr_list, rownames))
cat(sprintf("  三数据集基因交集: %d\n", length(common_genes)))

# 合并所有样本
all_expr <- do.call(cbind, lapply(gene_expr_list, function(x) x[common_genes, , drop = FALSE]))
all_group <- factor(unlist(group_list), levels = c("NN", "PP"))
all_batch <- factor(unlist(batch_list))

cat(sprintf("  合并矩阵: %d 基因 × %d 样本\n", nrow(all_expr), ncol(all_expr)))
cat(sprintf("  批次分布: %s\n", paste(levels(all_batch), table(all_batch), sep="=", collapse=", ")))
cat(sprintf("  分组分布: NN=%d, PP=%d\n", sum(all_group=="NN"), sum(all_group=="PP")))

# =============================================================================
# Step 5: ComBat 批次校正 (带生物学协变量 mod)
# =============================================================================
cat("\n=== Step 5: ComBat 批次校正 (mod 含 group) ===\n")

# mod 包含生物学协变量 (group)，防止 ComBat 把生物学差异当批次效应消除
mod <- model.matrix(~ all_group)
cat(sprintf("  mod 矩阵: %d 行 × %d 列\n", nrow(mod), ncol(mod)))

all_expr_combat <- ComBat(dat = all_expr, batch = all_batch, mod = mod)
cat("  ComBat 校正完成\n")

# =============================================================================
# Step 6: limma 差异分析 — eBayes(trend=TRUE)
# =============================================================================
cat("\n=== Step 6: limma eBayes(trend=TRUE) ===\n")

design <- model.matrix(~ all_group)
colnames(design) <- c("Intercept", "PP_vs_NN")

fit <- lmFit(all_expr_combat, design)
fit <- eBayes(fit, trend = TRUE)

deg_results <- topTable(fit, coef = "PP_vs_NN", number = Inf, adjust.method = "BH",
                        sort.by = "P")
deg_results$Gene <- rownames(deg_results)

# 添加调控方向
deg_results$regulation <- ifelse(deg_results$logFC > 0, "Up", "Down")

# =============================================================================
# Step 7: 阈值筛选
# =============================================================================
cat("\n=== Step 7: 阈值筛选 (|logFC|>1, FDR<0.05) ===\n")

deg_sig <- subset(deg_results, abs(logFC) > 1 & adj.P.Val < 0.05)
deg_sig <- deg_sig[order(deg_sig$adj.P.Val), ]

cat(sprintf("  DEG 总计: %d (Up=%d, Down=%d)\n",
            nrow(deg_sig), sum(deg_sig$regulation == "Up"),
            sum(deg_sig$regulation == "Down")))

# 检查之前缺失的基因是否进入 DEG
cat("  之前注释缺失的基因 DEG 状态:\n")
for (g in c("RUNX3", "ZNF683", "CAMP", "S100A7", "S100A8", "S100A9", "LCN2")) {
  idx <- which(deg_results$Gene == g)
  if (length(idx) == 0) {
    cat(sprintf("    %s: NOT in gene-level matrix\n", g))
  } else {
    dr <- deg_results[idx, ]
    in_deg <- abs(dr$logFC) > 1 & dr$adj.P.Val < 0.05
    cat(sprintf("    %s: logFC=%.3f, FDR=%.2e %s\n", g, dr$logFC, dr$adj.P.Val,
                ifelse(in_deg, "✅ DEG", "")))
  }
}

# =============================================================================
# Step 8: 输出结果
# =============================================================================
cat("\n=== Step 8: 保存结果 ===\n")

# 完整 DEG 结果 (所有基因)
write.table(deg_results,
            file = file.path(OUT_DIR, "results", "deg_results_all_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  全基因结果: deg_results_all_genes.tsv (%d 行)\n", nrow(deg_results)))

# 显著 DEG
dir.create(file.path(OUT_DIR, "results"), showWarnings = FALSE, recursive = TRUE)
write.table(deg_sig,
            file = file.path(OUT_DIR, "results", "deg_results_significant.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  显著 DEG: deg_results_significant.tsv (%d 行)\n", nrow(deg_sig)))

# ComBat 校正后表达矩阵
write.csv(all_expr_combat,
          file = file.path(OUT_DIR, "results", "expression_combat_corrected.csv"),
          quote = FALSE)
cat(sprintf("  表达矩阵: expression_combat_corrected.csv (%d 基因 × %d 样本)\n",
            nrow(all_expr_combat), ncol(all_expr_combat)))

# 样本信息
sample_info <- data.frame(
  sample = colnames(all_expr),
  group  = as.character(all_group),
  batch  = as.character(all_batch),
  stringsAsFactors = FALSE
)
write.table(sample_info,
            file = file.path(OUT_DIR, "results", "sample_info.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 注释映射表 (AnnoProbe)
write.table(probe2gene,
            file = file.path(OUT_DIR, "results", "probe2gene_AnnoProbe.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 运行摘要
sink(file.path(OUT_DIR, "results", "pipeline_summary_v2_AnnoProbe.txt"))
cat("=== 标准化 DEG 分析管线摘要 (v2 AnnoProbe) ===\n\n")
cat("数据集:\n")
for (ds in DATASETS) {
  cat(sprintf("  %s: %d 样本 (NN=%d, PP=%d)\n", ds,
              ncol(expr_list[[ds]]),
              sum(group_list[[ds]] == "NN"),
              sum(group_list[[ds]] == "PP")))
}
cat(sprintf("\n总样本: %d (NN=%d, PP=%d)\n",
            length(all_group), sum(all_group=="NN"), sum(all_group=="PP")))
cat(sprintf("基因交集: %d\n", length(common_genes)))
cat(sprintf("\nDEG (|logFC|>1, FDR<0.05): %d (Up=%d, Down=%d)\n",
            nrow(deg_sig), sum(deg_sig$regulation=="Up"), sum(deg_sig$regulation=="Down")))
cat(sprintf("FDR 范围: %.2e – %.2e\n", min(deg_sig$adj.P.Val), max(deg_sig$adj.P.Val)))
cat(sprintf("|logFC| 范围: %.2f – %.2f\n", min(abs(deg_sig$logFC)), max(abs(deg_sig$logFC))))
cat("\nTop 20 DEG:\n")
print(head(deg_sig[, c("Gene", "logFC", "AveExpr", "adj.P.Val", "regulation")], 20))
cat("\n方法引用:\n")
cat("  归一化: RMA (Irizarry et al. 2003, Biostatistics)\n")
cat("  批次校正: ComBat with mod (Johnson et al. 2007, Biostatistics)\n")
cat("  注释: AnnoProbe v0.1.0 idmap(type='bioc') (Zeng JM et al., Bioconductor hgu133plus2.db)\n")
cat("  多探针压缩: MaxMean (Miller et al. 2011, BMC Bioinformatics)\n")
cat("  差异检验: limma eBayes(trend=TRUE) (Ritchie et al. 2015, NAR)\n")
cat("  多重校正: Benjamini-Hochberg FDR (Benjamini & Hochberg 1995)\n")
sink()

cat("\n=== 管线完成 ===\n")
cat(sprintf("输出目录: %s/results\n", OUT_DIR))
