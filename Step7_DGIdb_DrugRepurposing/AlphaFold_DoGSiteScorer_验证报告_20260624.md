# SNP 含义 & AlphaFold/DoGSiteScorer 可用性验证报告

> 2026-06-24 | 响应问询 @19:23 GMT+8

---

## 问题1：文献结果表格中的 SNP 指什么？

**SNP = Single Nucleotide Polymorphism（单核苷酸多态性）**

在 Alghubayshi 2025（Front. Bioinform.）的上下文中：
- 基因来源是**沙特人群 SCD 的 GWAS**（Alghubayshi et al. 2025 前置工作）
- GWAS 芯片检测数十万 SNP 位点，显著位点（P < 5×10⁻⁸）→ 基因注释 → **31 个疾病严重程度相关基因**
- 表格中 SNP 列（rs编号如 rs2814778）即每个基因对应的 GWAS 显著变异位点

| 对比维度 | Alghubayshi 2025 | 我们的研究 |
|---------|-----------------|-----------|
| 基因来源 | GWAS（遗传关联） | DEG转录组 + PPARgene 功能交集 |
| 输入数量 | 31 基因 | 18 基因 |
| 问题类型 | "哪些遗传变异增加SCD严重性" | "哪些PPAR靶基因在银屑病中差异表达" |
| 可药性方法 | DGIdb v5 + DoGSiteScorer + AlphaFold | DGIdb v5 + 分子对接 |

两种策略回答不同层面问题，都合法。

---

## 问题2：AlphaFold + DoGSiteScorer 数据能否获取？

### ✅ AlphaFold DB — 完全可用

| Hub 基因 | UniProt ID | AlphaFold PDB | 文件大小 | 状态 |
|----------|-----------|--------------|---------|------|
| ADIPOQ | Q15848 | [AF-Q15848-F1-model_v6.pdb](https://alphafold.ebi.ac.uk/files/AF-Q15848-F1-model_v6.pdb) | <1 MB | ✅ |
| APOE | P02649 | [AF-P02649-F1-model_v6.pdb](https://alphafold.ebi.ac.uk/files/AF-P02649-F1-model_v6.pdb) | ~0.20 MB | ✅ |
| PLIN1 | O60240 | [AF-O60240-F1-model_v6.pdb](https://alphafold.ebi.ac.uk/files/AF-O60240-F1-model_v6.pdb) | <1 MB | ✅ |
| BCL2 | P10415 | [AF-P10415-F1-model_v6.pdb](https://alphafold.ebi.ac.uk/files/AF-P10415-F1-model_v6.pdb) | <1 MB | ✅ |

- 许可证：CC-BY 4.0（可自由学术使用）
- API：`https://alphafold.ebi.ac.uk/api/prediction/{UniProt_ID}`
- 所有结构均为 v6 版本（AlphaFold-Multimer 预测）

### ✅ DoGSiteScorer (ProteinsPlus) — 完全可用

**REST API 端点：** `https://proteins.plus/api/dogsite_rest`

**工作流程（已验证通过 APOE）：**
1. 上传 AlphaFold PDB 到 ProteinsPlus → 获得 custom PDB ID
2. 提交 DoGSite 任务（bindingSitePredictionGranularity=1 可获取 druggability 评分）
3. 轮询结果 → 获取每个口袋的 `drugScore` 和 `simpleScore`

**实际测试结果（APOE）：**
- 检测到 **41 个口袋**（含主口袋和子口袋）
- 顶层口袋 P_0：drugScore = **0.808**（>0.5 阈值 ✅）
- 子口袋分数范围：0.066 ~ 0.738

**限制：**
- 速率限制：30 任务/分钟
- DoGSiteScorer 有额外 CPU/RAM 使用限制（需分时段提交）
- 结果的 `drugScore` 在 0-1 范围，Alghubayshi 使用 >0.5 为可药性阈值

---

## 对我们论文的意义

### 可直接加入 §2.2.10 分子对接部分：

```
"为补充分子对接分析,我们使用 DoGSiteScorer (Volkamer et al., 2012, 
ProteinsPlus) 评估了 4 Hub 基因产物的结合口袋可药性。蛋白质结构
从 AlphaFold DB (Jumper et al., 2021) 获取，drugScore > 0.5 的
口袋被认为具有潜在可药性。"
```

### 可引用文献：
- Volkamer et al. 2012 (DoGSiteScorer)
- Jumper et al. 2021 (AlphaFold)
- Alghubayshi et al. 2025 (方法学参考)
- Callaway 2024 (2024 诺贝尔化学奖)

### ⚠️ 需注意：
- 我们的分子对接已经是对靶点结合能力的直接验证
- DoGSiteScorer 提供的是"口袋适合度"的物理化学评估，与对接互补但不重复
- 如果论文篇幅有限，DoGSiteScorer 可作为**补充方法**而非核心方法
- Alghubayshi 用它是因为他们研究"novel targets"（无已知药物的基因），我们需要它的情况不同（4个Hub中3个有已知药物）
