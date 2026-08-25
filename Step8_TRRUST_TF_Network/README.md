# Step8_TRRUST_TF_Network

## 目录结构

```
Step8_TRRUST_TF_Network/
├── README.md                                    # 本文件
├── raw_data/
│   └── trrust_rawdata_human.tsv                 # TRRUST v2 原始数据（人类，9,396条）
├── data/
│   ├── TRRUST_enrichment_significant_human.tsv  # 人类TF显著富集（12个，Table 12来源）
│   ├── trrust_interactions_aggregated.tsv       # TF→Hub聚合边（54条，含n_PMIDs）
│   ├── trrust_edges.tsv                         # 网络边列表
│   └── trrust_nodes.tsv                         # 网络节点列表（4 Hub + 49 TF）
├── figures/
│   ├── TRRUST_TF_network_4hubs.png              # TF-Hub网络图（PNG）
│   └── TRRUST_TF_network_4hubs.pdf              # TF-Hub网络图（PDF，论文用）
└── scripts/
    └── build_trrust_network.py                  # 构建4Hub TF网络脚本
```

## 数据说明

### TRRUST v2 数据库

- **来源**: Han H et al. Nucleic Acids Res. 2018;46(D1):D380-D386. PMID: 29087512
- **下载时间**: 2026-06-21
- **下载地址**: https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
- **版本**: TRRUST v2（最终版本，无v3）
- **人类TF-靶基因对**: 9,396条

### 富集分析（data/TRRUST_enrichment_significant_human.tsv）

- **方法**: Enrichr API（library=TRRUST_Transcription_Factors_2019），BH-FDR<0.05
- **输入**: 18个PPAR-targeted DEGs
- **人类显著TF**: 12个
- **Top 3**: PPARA (FDR=5.42e-04), PPARG (FDR=1.12e-03), SP1 (FDR=1.52e-03)

### TF→Hub基因网络（data/trrust_interactions_aggregated.tsv）

- **节点**: 53个（4 Hub + 49 TF）
- **聚合边**: 54条
- **Hub边数分布**: ADIPOQ=3, APOE=3, PLIN1=2, BCL2=46
- **关键TF**: NFKB1/RELA→PLIN1(Activation), PPARG→ADIPOQ+BCL2(Mixed), STAT3→BCL2(Mixed)

## 论文中的用途

| 用途 | 文件 |
|------|------|
| Table 12: Top人类显著TF | `data/TRRUST_enrichment_significant_human.tsv` |
| Figure 15: TF-Hub网络图 | `figures/TRRUST_TF_network_4hubs.pdf` |
| Results正文数据 | `data/trrust_interactions_aggregated.tsv` |
| Methods: TRRUST版本+参数 | `raw_data/trrust_rawdata_human.tsv` |

## 注意事项

1. **BCL2边数偏多**: 46/54条边连接BCL2（凋亡核心基因），论文中需说明此不平衡
2. **TRRUST v2时效性**: 数据库最后更新2018年，Methods中需注明版本
3. **Mixed/Unknown占比较高**: 因TRRUST整合不同实验条件结果，论文中需说明context-dependent特性
