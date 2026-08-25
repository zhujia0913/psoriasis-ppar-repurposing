# DGIdb 药物重定位管线 — 严谨方法论

## 参考文献基础

| 文献 | 期刊 | 方法论贡献 |
|------|------|-----------|
| **Freshour SL et al. 2021** (PMID: 33237278) | *Nucleic Acids Res* | DGIdb 4.0 架构：41数据源、药物标准化(ChEMBL+Wikidata)、Query/Interaction Score 算法、交互方向性框架 |
| **Chen Q et al. 2025** (PMID: 41511918) | *未标期刊* | DEG→PPI→MCODE→CytoHubba→TRRUST→DGIdb 标准管线，过敏性鼻炎药筛方法学基准 |
| **Li X et al. 2025** (PMID: 40406126) | *Front Immunol* | 银屑病-克罗恩病共享生物标志物+DEG→WGCNA→DGIdb 药筛，2025年发表 |
| **Fan J et al. 2024** (PMID: 39072329) | *Front Immunol* | 银屑病-动脉粥样硬化铁死亡/坏死性凋亡 DEG→WGCNA→候选药物预测 |
| **Griffith M et al. 2013** (PMID: 24122041) | *Nat Methods* | DGIdb 原始发表 |

---

## 方法学管线（8 步，4 层）

```
Layer 1: 靶点识别          Layer 2: 药物-靶点映射       Layer 3: 药物过滤           Layer 4: 临床验证
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│ Step 0: DEG鉴定  │     │                 │     │                 │     │                 │
│ Step 0a: Hub筛选 │ ──→ │ Step 1: DGIdb API│ ──→ │ Step 5: 药物类别  │ ──→ │ Step 7: FDA验证  │
│ Step 0b: 表达验证│     │ Step 2: DEG验证  │     │ Step 6: 安全排除  │     │ Step 8: ATC分类  │
│                  │     │ Step 3: 临床状态  │     │                 │     │                 │
│                  │     │ Step 4: 非药物    │     │                 │     │                 │
└─────────────────┘     └─────────────────┘     └─────────────────┘     └─────────────────┘
```

---

### Layer 1: 靶点识别（论文方法部分）

#### Step 0: 差异表达基因（DEG）鉴定

**方法**: GEO 数据集 → AnnoProbe 重注释 → ComBat 批次校正 → limma eBayes(trend=TRUE)

| 参数 | 值 |
|------|-----|
| 数据集 | GSE13355, GSE14905, GSE78097（209样本: 106 NN + 103 PP）|
| 平台 | GPL570 (Affymetrix HG-U133 Plus 2.0) |
| 注释方法 | AnnoProbe `type="pipe"`（20,188 基因池） |
| 批次校正 | ComBat（参考批次=GSE14905） |
| 显著性阈值 | FDR < 0.05 且 |logFC| > 1 |
| 软件 | R 4.5.0, GEOquery, sva, limma, AnnoProbe |

**关键决策**: 采用 AnnoProbe 而非 GPL570 官方注释表，基因池从 10,226 扩展到 20,188（+97.4%），DEG 数量从 352 增至 759。

#### Step 0a: Hub 基因筛选

**方法**: 疾病相关基因集（PPARgene v2.0 Human-only, 177 基因）∩ DEG（759）→ 18 交集基因 → STRING PPI 网络（confidence ≥ 0.4）→ CytoHubba 11 算法 → Rank Sum 截断

| 步骤 | 工具/参数 |
|------|-----------|
| 疾病基因集 | PPARgene v2.0 (Fang et al. 2016, *Nucleic Acids Res*), 仅保留物种=Human（177 基因） |
| 交集筛选 | PPARgene ∩ DEG = 18 基因 |
| PPI 网络 | STRING v12.0, confidence ≥ 0.4, 物种=Homo sapiens |
| Hub 筛选 | CytoHubba v0.1 (Chin et al. 2014, *BMC Syst Biol*), 11 算法全跑 |
| 截断方法 | Rank Sum（Degree + Betweenness + MNC + Closeness, 4 算法），自然断点截断 |
| 最终 Hub | ADIPOQ, APOE, PLIN1, BCL2（全部下调） |

**关键决策**: 弃用 ≥8/11 交集法（无文献依据）和完全交集法（导致空集），采用 Rank Sum + 自然断点（Li et al. 2020, *Front Genet* 支撑）。

#### Step 0b: 组织表达验证

**方法**: Human Protein Atlas (HPA) + GTEx + FANTOM5

| 验证维度 | 标准 | 排除条件 |
|---------|------|---------|
| 皮肤表达 | HPA/GTEx 皮肤组织 nTPM > 0 | 无皮肤表达 |
| 睾丸蛋白 | HPA 蛋白 IHC 生殖细胞阳性 | 生殖细胞阳性（生殖毒性风险） |
| 组织特异性 | GTEx 组织特异性评分 | 无此标准 |

**结果**: 4 Hub 基因全部通过（PLIN1 睾丸仅 Sertoli 细胞阳性，非生殖细胞）。

---

### Layer 2: 药物-靶点映射

#### Step 1: DGIdb API 查询 + Flag 分类

**方法**: DGIdb v5.0 GraphQL API → 按 concept_id 和药物名分类标记

**Concept ID 前缀过滤**（API返回153条总交互，分布如下）:

| concept_id 前缀 | 含义 | 条数 | 处理 |
|----------------|------|------|------|
| `rxcui:` | FDA批准药物 (RxNorm) | 67 | **保留** |
| `chembl:` | 研究化合物 | 4 | 排除 |
| `ncit:` | NCI Thesaurus 药物类别名/研究试剂 | 46 | 排除 |
| `iuphar.ligand:` | IUPHAR 配体（研究工具） | 31 | 排除 |
| `drugbank:` | DrugBank 条目 | 4 | 排除 |
| `wikidata:` | Wikidata 条目 | 1 | 排除 |

**仅保留 rxcui-tagged 交互**（67条，对应65个unique drugs；FENOFIBRATE MICRONIZED和TRETINOIN各出现2次，分别靶向不同Hub基因）。

后续 Flag 分类在 rxcui 范围内进行:

| Flag | 定义 | 处理 |
|------|------|------|
| OK | 可进入后续步骤 | 进入 Step 2 |
| CHEMBL | `concept_id` 以 `chembl:` 开头（研究化合物） | 排除 |
| RESEARCH_TOOL | IUPHAR 配体编号 + GtP/NCI 来源（研究工具） | 排除 |
| SUPPLEMENT | 维生素/矿物质/营养补充剂 | 排除 |
| DRUG_CLASS | NCIT 药物类别名（非特定药物） | 排除 |

**DGIdb API 调用参数**:
- 端点: GraphQL `https://dgidb.org/api/graphql` (DGIdb v5.0, 2026-06-24 查询)
- 返回字段: `gene`, `drug`, `concept_id`, `interaction_types`, `score`, `sources`, `n_pmids`, `pmids`
- 基因输入: ADIPOQ, APOE, PLIN1, BCL2

**Interaction Score 计算**（基于 Freshour et al. 2021, Figure 3）:
```
Evidence_Score = ln(n_pmids + 1) + ln(n_sources + 1)
Interaction_Score = Evidence_Score × (1 / ln(e + n_drugs_per_gene)) × (1 / ln(e + n_genes_per_drug))
```
惩罚高 promiscuity 基因/药物的非特异性交互。

**Directionality 标注**（基于 Freshour et al. 2021, Supp. Table S2）:
- agonist/activator/positive modulator → ACTIVATING
- inhibitor/antagonist/negative modulator → INHIBITORY
- 其他/unknown → N/A
- 对照 Hub 基因下调方向: ACTIVATING = MATCH, INHIBITORY = MISMATCH

#### Step 2: DEG 表达一致性验证

**方法**: 所有 Hub 基因已在 DEG 列表中（必要条件），按 logFC 方向和显著性交叉验证。

| 验证项 | 标准 |
|--------|------|
| 靶基因是否在 DEG 列表 | 必须 ≥ 1 个 Hub 基因是 DEG |
| 方向一致性 | 记录每个药物-基因对的 logFC 值 |
| HPA/GTEx 交叉验证 | 组织表达模式与疾病组织匹配 |

#### Step 3: 临床状态排除

**方法**: 从 FDA Orange Book、Drugs@FDA、EMA 撤回列表中排除已撤回/终止的药物。

| 排除类别 | 示例 |
|---------|------|
| Withdrawn（安全性） | Troglitazone (肝毒性), Rofecoxib (心血管) |
| Terminated（无效） | Torcetrapib (Phase 3 死亡率增加) |
| Safety-restricted | Rosiglitazone (FDA 黑框警告) |

#### Step 4: 非药物产品排除

**方法**: 排除研究工具肽、内源性代谢物、非药物化学品。

| 排除类别 | 示例 |
|---------|------|
| PAR2 研究工具 | AC-55541, SLIGKV-NH2 |
| 内源性分子 | Heme, C5A, Lysozyme |
| 非药物化学品 | Citric Acid, Hydroquinone, Staurosporine |
| 抗体-药物偶联物 | Brentuximab Vedotin（复杂生物制剂）|

---

### Layer 3: 药物过滤

#### Step 5: 药物类别排除

**原则**: 排除临床不可行或将引入显著混杂偏倚的类别。

| 排除类别 | 说明 | 示例 |
|---------|------|------|
| CHEMOTHERAPY | 抗肿瘤化疗药，重度毒性不可接受 | Paclitaxel, Doxorubicin, Venetoclax |
| ALZHEIMERS | 适应症不匹配，CNS 特异性 | Donepezil, Galantamine |
| HORMONE_CONTRACEPTIVE | 激素类药物，干扰研究 | Norethindrone Acetate |
| ANTIVIRAL/ANTIMICROBIAL | 抗感染药，适应症不匹配 | Nevirapine, Ganciclovir |
| NON_DRUG | 非药物化学品/代谢物 | Beauvericin (真菌毒素) |
| PERIPHERAL | 外周作用药物，无法系统重定位 | Carbetapentane (止咳药) |
| OBSOLETE | 已被更安全药物替代 | Clofibrate → Fenofibrate |
| INVESTIGATIONAL | 早期临床阶段，安全性数据不足 | 4-Phenylbutyric Acid |

#### Step 6: 安全性排除

**原则**: 基于已知副作用数据库和银屑病加重证据。

| 排除类别 | 证据来源 | 示例 |
|---------|---------|------|
| ANTICONVULSANT | FDA 标签 - CNS 副作用 | Phenytoin, Carbamazepine |
| ANTIDEPRESSANT | FDA 标签 - 血清素综合征风险 | Fluoxetine, Sertraline |
| GLUCOCORTICOID | 长期免疫抑制 + 银屑病反弹 | Prednisone, Dexamethasone |
| OPIOID | 成瘾性 + FDA 管制 | Morphine, Fentanyl |
| PSORIASIS_AGGRAVATOR | PubMed 病例报告/综述 (PMID: 19389576 等) | Lithium, Hydroxychloroquine, Propranolol, ACE 抑制剂 |
| ANTITHYROID | 非适应症，甲状腺毒性 | Methimazole |

---

### Layer 4: 临床验证

#### Step 7: FDA 审批状态

**方法**: Drugs@FDA / Orange Book 查询审批状态与申请号。

| 状态 | 标注 |
|------|------|
| NDA | 新药申请号 |
| ANDA | 仿制药申请号 |
| BLA | 生物制剂许可申请 |
| OTC | 非处方药 |
| Pending | 未确认 |

#### Step 8: WHO ATC 分类

**方法**: WHO Collaborating Centre for Drug Statistics Methodology ATC/DDD Index。

| ATC 前缀 | 治疗类别 |
|---------|---------|
| C10 | 脂质调节剂 |
| D05 | 抗银屑病药 |
| D10 | 抗痤疮药 |
| A10 | 糖尿病用药 |
| C03 | 利尿剂 |
| C09 | RAS 抑制剂 |
| A05 | 胆汁酸制剂 |

---

### 候选药排序标准

最终候选药按 DGIdb Interaction Score 降序排列，无分层分级。

**Interaction Score** = Evidence Score × Relative Drug Specificity × Relative Gene Specificity（DGIdb v5 定义，详见 GitHub Issue #515）。Score 越高表示药物-基因交互证据越强、特异性越高。

---

### 产出文件清单

| 文件 | 内容 |
|------|------|
| `step1_dgidb_raw.tsv` | DGIdb GraphQL 原始查询结果（153 条交互）|
| `step1_6_full.tsv` | 通过 6 步过滤的全部 86 个药物（含过滤标签、IScore、Directionality，20 列）|
| `final_candidates.tsv` | 最终 11 候选药交互（按 Interaction Score 降序排列）|
| `4hub_hpa_expression_summary.tsv` | 4 Hub 基因 HPA/GTEx 组织表达汇总 |

**查询方式**: DGIdb v5.0 网页 / GraphQL API 查询（查询脚本未随仓库分发；原始响应保留于本文件夹 `dgidb_v5_4hub_raw_response.json`，交互明细见 `dgidb_v5_4hub_all_interactions.tsv`，最终候选见 `DGIdb_v5_final_candidates.tsv`）
- R3 验证: 每步输出 `f-string` 变量化（禁止硬编码数字）
- 交叉验证: 最终输出 logFC vs `PPI_nodes.tsv` 全量 diff

---

### 方法学局限与透明性声明

1. **DGIdb interaction_type 覆盖率**: ~76% 交互来自 NCI/PharmGKB 源，返回 `interaction_type=none`。DGIdb 4.0 论文确认此限制（仅 ~24% 交互有 curated directionality）。本研究中候选药全部来自 NCI/PharmGKB 源，因此 Directionality 全部为 N/A，不影响最终排序。

2. **BCL2 药物低置信度**: Isotretinoin, Tazarotene, Tretinoin, Ursodiol 的 BCL2 关联来自 NCI therapeutic classification，DGIdb Score < 0.5 且无具体 interaction type，标注为 LOW 置信度（Discussion 中注明）。

3. **Interaction Score 动态性**: Freshour et al. 论文指出随数据源更新，Interaction Score 可变。

4. **网络邻近性未纳入**: 因 15 节点网络下网络邻近性分析统计效力不足（需 50+ 节点网络），本管线未包含 Z-score 验证，作为局限性标注。
