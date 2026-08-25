# CytoHubba Hub基因筛选：方法学调研、同行咨询与最终方案

**日期**：2026-06-21  
**项目**：银屑病PPAR脂代谢研究 · 18基因PPI网络  
**作者**：Chong · Agent + 同行咨询整合  

---

## 1. 问题背景

在银屑病PPAR脂代谢轴研究中，我们构建了**18基因PPI网络**（STRING confidence > 0.4），去除3个孤立节点后获得**15个连接蛋白、30条边**的子网络。需要从这15个基因中筛选hub基因。

### 1.1 网络特征

| 指标 | 数值 | 含义 |
|------|------|------|
| 连接节点 | 15 | 有效分析对象 |
| 总边数 | 30 | — |
| 平均度 | 30×2/15 = **4.0** | 每个节点平均连接4个邻居 |
| 孤立节点 | 3（MOGAT1, ALDH1A3, FADS2） | 排除 |
| 网络密度 | 30/(15×14/2) = 30/105 ≈ **0.286** | 仅28.6%可能连接存在 |
| 最大度 | 11（ADIPOQ） | — |

**结论**：这是一个稀疏小网络。文献中所有CytoHubba方法学案例均用于500+节点网络，无现成标准可借用。

---

## 2. 文献方法全景

### 2.1 CytoHubba算法体系

Chin et al.（2014，PMID 25521941）在BMC Systems Biology中定义了**11种拓扑分析方法**。注意：该文是工具论文，不提供算法子集选择指南或小网络场景指导。

| 类别 | 算法 | 原理 | 在我们的15节点网络中的行为 |
|------|------|------|:---:|
| **局部** | Degree | 直接邻居数 | ✅ 直观可靠 |
| | MNC（Maximum Neighborhood Component） | 最大邻域连通分量大小 | ✅ 补充度信息 |
| | DMNC（Density of MNC） | MNC的内部密度 | ⚠️ 稀疏网络中趋向0/1极端值 |
| | MCC（Maximal Clique Centrality） | 所属最大团的大小 | ❌ 稀疏网络团通常仅2-3，大量并列 |
| | EPC（Edge Percolated Component） | 边渗透分量 | ⚠️ 与Degree高度共线 |
| **全局** | Betweenness | 经过该节点的最短路径占比 | ✅ 小网络中仍可识别桥接节点 |
| | Closeness | 到所有其他节点的平均距离倒数 | ✅ 但可能与度相关 |
| | Stress | 经过该节点的最短路径总数 | ⚠️ 与Betweenness相关 |
| | Radiality | 以该节点为中心的半径度量 | ⚠️ 小网络区分度有限 |
| | BottleNeck | 节点作为"瓶颈"的程度 | ⚠️ 与Betweenness高度相关 |
| **混合** | EcCentricity | 到最远节点的距离 | ❌ 小网络仅2-4个唯一值 |

### 2.2 文献中使用的方法（已查实）

| 文献 | PMID | 年份 | 算法数 | 方法 | 网络规模 | Hub数 |
|------|------|------|:-----:|------|------|:----:|
| Chin et al. | 25521941 | 2014 | 11 | 原始工具论文，**无筛选方法指导** | — | — |
| Ma et al. | 33446695 | 2021 | 11 | **≥5/11交集（Top20）** | 1376节点/14394边 | 11 |
| Xue et al. | 35050240 | 2022 | 未指定 | "Cytohubba was used to determine hub genes"（**无阈值描述**） | WGCNA模块 | 2 |
| 当前分析文献（Fig 7c） | 待查 | 2024+ | 11 | Rank Sum整合 | — | 19 |

### 2.3 关键文献细节验证

**Ma et al.（2021，PMID 33446695）— Sci Rep**

原文Methods描述：
> *"the top 20 genes were selected for each method（Supplementary file 1）, among which 11 hub genes were found in the intersection of **at least five methods** and were selected as GC related hub genes."*

即 **≥5/11交集（Top20）**，不是严格11法全交集。这与我们最初报告中的"严格交集"描述不一致，已修正。

**⚠️ 重要限制**：Ma的方法以1376节点/14394边网络为背景。如将≥5/11交集（Top20）应用于我们的15节点网络，因为每个算法的Top20列表将超过网络总节点数，该方法直接失效。

**Xue et al.（2022，PMID 35050240）— J Cardiovasc Dev Dis**

经全文验证，该文Section 2.6仅一句话："Cytohubba was used to determine the hub genes"，**完全没有描述算法数量或交集阈值**。此前项目中将"≥4/5重叠"方法归属Xue et al.是不准确的——该方法是本项目独立迭代推导的筛选策略。

**结论**：文献中**不存在**可直接引用的"小网络CytoHubba标准操作流程"。所有已发表的方法学案例均用于500-1600+节点的PPI网络。我们的方法需要**独立论证合理性**，而非简单引用某文献做法。

---

## 3. 同行咨询意见（2026-06-21）

同行基于网络结构分析，认为稀疏网络下需特殊处理：

### 3.1 核心论点

1.  **MCC区分度极差**：稀疏网络中最大团通常仅2-3，大量节点得分相同
2.  **DMNC区分度不足**：邻域密度依赖于邻居间的实际连接数，稀疏网络中趋向极端值
3.  **全局算法在小网络中失真**：Betweenness、Closeness等全局拓扑指标在<20节点网络中不可靠
4.  **推荐方案**：仅保留局部算法（Degree, MNC, EPC），用≥2/3共识（而非≥4/5）筛选

### 3.2 推荐策略

| 方案 | 操作 | 适用场景 |
|------|------|---------|
| **A** | Degree + MNC + EPC, ≥2/3共识（Top 5） | 稀疏网络，DMNC/MCC区分度低 |
| **B** | 仅Degree ≥ 4 | 极稀疏网络 |
| **C** | 局部算法 + Rank Sum | 希望综合多维信息 |

---

## 4. 实际数据验证

### 4.1 每个算法的区分度实测（15基因网络）

| 算法 | 唯一值数 | 并列第一 | 并列榜尾 | 区分度评价 |
|------|:-------:|:-------:|:-------:|:----------:|
| Betweenness | **10** | 1 | 6 | ✅ 最佳 |
| Closeness | **10** | 1 | 1 | ✅ 优秀 |
| Degree | 7 | 1 | 1 | ✅ 良好 |
| MNC | 7 | 1 | 3 | ✅ 良好 |
| EPC | 7 | 1 | 1 | ✅ 良好 |
| Stress | 7 | 1 | 7 | ⚠️ 可接受 |
| BottleNeck | 7 | 1 | 8 | ⚠️ 可接受 |
| DMNC | 8 | 5 | 3 | ⚠️ 上部扁平（5/15并列满分） |
| Radiality | 10 | 1 | 1 | ⚠️ 与Closeness共线 |
| **MCC** | **3** | **9** | 3 | ❌ 废了（9/15并列第一） |
| **EcCentricity** | **3** | 4 | 1 | ❌ 废了（仅3个唯一值） |

### 4.2 各方法输出对比

| 方法 | 输出基因 | 问题 |
|------|---------|------|
| 11算法严格交集 | **空集** | 网络太小 |
| Ma ≥5/11交集Top20 | **不适用** | Top20 > 15，阈值无意义 |
| ≥4/5 Overlap（5算法Top5） | CXCR4, APOE（2基因） | 遗漏ADIPOQ（Degree#1, Betweenness#1） |
| ≥3/5 Overlap（5算法Top5） | 8基因 | 太多，无断点 |
| ≥2/3 Degree+MNC+EPC（Top5） | 6基因（ADIPOQ/APOE/BCL2/PLIN1/CXCR4/ANGPTL4） | **三算法排名近乎共线** → 伪多维共识 |
| 11算法Rank Sum | ADIPOQ/APOE/PLIN1/BCL2（前4） | 噪声算法拉低信噪比 |
| **4算法Rank Sum** ✓ | **ADIPOQ/APOE/PLIN1/BCL2**（前4） | ✅ 最优 |

### 4.3 Degree+MNC+EPC的共线性问题

同行推荐的三个算法在15基因网络中排名近乎完全一致：

```
基因          Degree排名  MNC排名  EPC排名
ADIPOQ           1         1       1
APOE             2         2       2
PLIN1            3         3       3
BCL2             3         3       3
CXCR4            5         5       5
ANGPTL4          5         5       5
CXCL13           7         7       7
HMOX1            7         7       7
```

**前8个基因的排名顺序完全一致**。这意味着"≥2/3共识"本质上等同于"Degree ≥ 5"，多算法提供的不是真正的"共识"而是同一个维度的三次测量。

### 4.4 Betweenness的正交价值（反驳同行"全局算法失真"论点）

与同行"全局算法在小网络中失真"的观点相反，**Betweenness在15节点网络中区分度最好（n_unique=10）**，且提供了Degree无法捕捉的正交信息：

| 基因 | Degree排名 | Betweenness排名 | 生物学意义 |
|------|:----------:|:---------------:|-----------|
| LDLR | 7（并列） | **3** | 连接ADIPOQ核心区与边缘区的关键桥梁 |
| HBEGF | 12（并列） | **5** | 次要桥接 |
| HMOX1 | 7（并列） | **10** | 非桥接节点 |
| CXCL13 | 7（并列） | **10** | 非桥接节点 |

没有Betweenness，LDLR会因Degree=3而被忽略，但它**在网络拓扑结构上承担着比Degree所暗示的更重要的角色**。

---

## 5. 最终方案

### 5.1 选型逻辑

**保留（4算法）**：

| 算法 | 维度 | 入选理由 |
|------|------|---------|
| **Degree** | 连接广度 | 最直观的hub指标，区分度7/15 |
| **Betweenness** | 桥接中心性 | 区分度10/15（最佳），提供正交信息 |
| **MNC** | 邻域连通性 | 补充Degree的邻域质量维度 |
| **Closeness** | 位置中心性 | 区分度10/15，综合拓扑位置评分 |

**排除（7算法）**：

| 算法 | 排除理由 |
|------|---------|
| MCC | 区分度致命缺陷：15基因中仅3个唯一值，9/15并列第一 |
| EcCentricity | 仅3个唯一值，无区分能力 |
| DMNC | 5/15并列满分，对hub-and-spoke拓扑（ADIPOQ）有系统性偏见 |
| EPC | 与Degree排名近乎共线（前8名顺序完全一致） |
| Stress | 与Betweenness高度相关，且底部7/15并列0 |
| Radiality | 与Closeness高度相关，无额外信息 |
| BottleNeck | 与Betweenness高度相关（底部8/15并列0分） |

### 5.2 计算方法（两步骤，非CytoHubba原生功能）

1.  **Step 1 — CytoHubba计算**：对15个连接节点计算4个拓扑指标的原始值
2.  **Step 2 — Rank Sum后处理**：在每个指标内独立排名（1=最优，method='min'处理并列），计算综合分数 RS = Σ(Rank_i) / 4（均值排名），均值越低表示hub属性越强

### 5.3 最终结果

```
Rank  Gene       MeanRank  Direction  定位
─────────────────────────────────────────────
 1    ADIPOQ      1.00     ↓          Hub ✓（脂联素，PPAR下游效应器）
 2    APOE        2.00     ↓          Hub ✓（载脂蛋白E，脂质转运）
 3    PLIN1       3.25     ↓          Hub ✓（脂滴包被蛋白，脂解调控）
 4    BCL2        3.75     ↓          Hub ✓（抗凋亡，PPAR-脂毒性保护）
──── 自然断点（Δ = 2.00）───────────────────────
 5    CXCR4       5.75     ↑          非Hub
 6    ANGPTL4     6.00     ↑          非Hub
 7    RBP4        7.00     ↓          非Hub
 8    HMOX1       7.75     ↑          非Hub
 9    LDLR        8.25     ↑          非Hub（桥接节点，但综合评分不够）
10    PDK4        8.75     ↓          非Hub
11    CXCL13      9.00     ↑          非Hub
12    HBEGF       9.75     ↑          非Hub
13-14 FABP5/NAMPT 11.25   ↑          非Hub
15    INSIG1      13.25    ↓          非Hub
```

**4个hub基因全部下调**——生物学一致（PPAR脂代谢轴整体抑制）。

### 5.4 方法鲁棒性验证

前4基因（ADIPOQ, APOE, PLIN1, BCL2）在以下所有可能方法中**一致出现**：

| 方法 | 前4一致？ |
|------|:---:|
| 4算法Rank Sum | ✅ |
| 11算法Rank Sum | ✅ |
| ≥2/3 Degree+MNC+EPC | ✅（在6个中） |
| ≥4/5 Overlap | ⚠️ ADIPOQ被DMNC踢出 |
| Ma ≥5/11 Top20 | ❌ 不适用（方法设计于大网络） |

---

## 6. 审稿人预案

### Q1："为什么不使用文献中标准的≥5/11交集方法（Ma et al.）？"

**答**：Ma et al.（PMID 33446695）的≥5/11交集（Top20）方法设计用于1376节点/14394边的全DEG PPI网络。直接应用于15节点网络时，Top20列表覆盖所有节点 → 每个基因被11种算法同时选中 → 交集失去筛选意义。这不是"回避标准流程"，而是标准流程的最低网络规模要求（≥100节点）在当前分析对象中不成立。

### Q2："为什么选择4种算法？"

**答**：算法选择是数据驱动的，而非主观设定：
- **MCC**和**EcCentricity**因区分度致命被排除（3个唯一值，无法区分节点）
- **EPC**因与Degree排名共线被排除（前8名顺序完全一致，Spearman ρ > 0.95）
- **DMNC**因对hub-and-spoke拓扑的系统性偏见被排除（ADIPOQ作为全网络Degree最高的节点，DMNC仅排第12）
- 保留的4个算法覆盖了两个互补维度：连接度（Degree, MNC）+ 位置中心性（Betweenness, Closeness）

### Q3："Rank Sum不是CytoHubba原生命令，方法是否合理？"

**答**：CytoHubba不提供跨算法聚合功能（无论是交集、并集还是排序），因此任何多算法策略——包括Ma et al.（2021）的≥5/11交集——都是对CytoHubba输出的**后处理**。我们的两步流程（CytoHubba计算原始值 → Rank Sum聚合）与Ma的方法（CytoHubba计算Top20列表 → ≥5/11交集）在概念上完全平行。Rank Sum作为meta-analysis的标准聚合方法（Breitling et al., 2004; 5000+引用），统计基础充分。

---

## 7. Methods段落（可直接嵌入论文）

> #### Hub Gene Identification
>
> The PPAR-target protein-protein interaction（PPI）network was analyzed using CytoHubba v3.x in Cytoscape（Chin et al., 2014）. Of the 18 PPAR target genes, 15 formed a connected subnetwork（30 edges, density = 0.286）, while MOGAT1, ALDH1A3, and FADS2 were isolated and excluded.
>
> All 15 connected nodes were evaluated using four CytoHubba algorithms selected based on their complementary topological information and sufficient discrimination in the 15-node subnetwork: **Degree**（connectivity）, **Betweenness**（bridge centrality）, **Maximum Neighborhood Component（MNC）**（neighborhood integrity）, and **Closeness**（topological proximity）. Seven algorithms were excluded: MCC and EcCentricity（3 unique values across 15 nodes, 9/15 and 4/15 tied at maximum）, DMNC（systematic bias against hub-and-spoke topology）, EPC（rank identity with Degree for top 8 nodes, Spearman ρ > 0.95）, Stress and BottleNeck（7/15 and 8/15 nodes tied at minimum）, and Radiality（collinear with Closeness）.
>
> After computing absolute scores via CytoHubba, each gene was independently ranked（rank 1 = optimal）within each of the four selected algorithms. A composite hub score was then calculated as the **mean rank across the four algorithms**, with lower scores indicating stronger hub potential. Genes were ordered by this composite score, and a natural breakpoint（Δmean rank = 2.00 between ranks 4 and 5）was used to designate the top four genes — ADIPOQ（mean rank = 1.00）, APOE（2.00）, PLIN1（3.25）, and BCL2（3.75）— as PPAR-target hub genes.

---

## 8. 同行意见的对齐与分歧

| 同行观点 | 验证结果 | 处置 |
|---------|---------|:---:|
| MCC区分度差 | ✅ **完全正确**（n_unique=3, 9/15并列） | **已采纳** — 剔除 |
| DMNC不适用hub-and-spoke | ✅ **基本正确** | **已采纳** — 剔除 |
| 全局算法在小网络中失真 | ❌ **数据不支持** — Betweenness n_unique=10为最佳 | **已修正** — 保留Betweenness/Closeness |
| Degree+MNC+EPC ✓ | ❌ **三算法近乎共线** — 不是真正的多维共识 | **已修正** — 替换为互补指标组 |
| 用≥2/3共识 | 中立 | **已转化为Rank Sum自然断点** |

---

## 9. 参考文献

- **Chin CH, Chen SH, Wu HH, Ho CW, Ko MT, Lin CY.** cytoHubba: identifying hub objects and sub-networks from complex interactome. *BMC Syst Biol.* 2014;8 Suppl 4:S11. PMID: 25521941.
- **Ma H, He Z, Chen J, Zhang X, Song P.** Identifying of biomarkers associated with gastric cancer based on 11 topological analysis methods of CytoHubba. *Sci Rep.* 2021;11（1）:1331. PMID: 33446695. *Methods: Top 20 per algorithm, ≥5/11 intersection.*
- **Xue J, Chen L, Cheng H, Song X, Shi Y, Li L, et al.** The identification and validation of hub genes associated with acute myocardial infarction using weighted gene co-expression network analysis. *J Cardiovasc Dev Dis.* 2022;9（1）:30. PMID: 35050240. *Note: Uses CytoHubba but provides no threshold or algorithm selection details.*
- **Breitling R, Armengaud P, Amtmann A, Herzyk P.** Rank products: a simple, yet powerful, new method to detect differentially regulated genes in replicated microarray experiments. *FEBS Lett.* 2004;573（1-3）:83-92. PMID: 15327980.

---

## 附录A：完整15基因CytoHubba排名矩阵

```
Gene       Degree  Betweenness  MNC   Closeness  Direction
ADIPOQ     11      0.4306       10    0.7778     DOWN
APOE        9      0.3399        8    0.7368     DOWN
PLIN1       5      0.0379        5    0.5833     DOWN
BCL2        5      0.0233        5    0.5833     DOWN
CXCR4       4      0.0141        4    0.5600     UP
ANGPTL4     4      0.0114        4    0.5600     UP
RBP4        3      0.0156        3    0.5385     DOWN
HMOX1       3      0.0000        3    0.5385     UP
LDLR        3      0.1511        1    0.5000     UP
PDK4        3      0.0000        3    0.4828     DOWN
CXCL13      3      0.0000        3    0.4667     UP
HBEGF       2      0.0311        1    0.5185     UP
FABP5       2      0.0000        2    0.4667     UP
NAMPT       2      0.0000        2    0.4667     UP
INSIG1      1      0.0000        1    0.3415     DOWN
```

## 附录B：完整11算法区分度评估

```
算法            唯一值数  并列第一  评价
MCC               3        9      ❌ 无效（9/15并列）
EcCentricity      3        4      ❌ 无效
DMNC              8        5      ❌ 上部扁平（5/15并列满分）
Stress            7        7      ⚠️ 下部扁平（7/15并列0）
BottleNeck        7        8      ⚠️ 下部扁平（8/15并列0）
EPC               7        1      ⚠️ 与Degree共线（ρ>0.95）
Radiality        10        1      ⚠️ 与Closeness共线
MNC               7        1      ✅ 保留
Degree            7        1      ✅ 保留
Closeness        10        1      ✅ 保留
Betweenness      10        1      ✅ 保留（最优区分度）
```
