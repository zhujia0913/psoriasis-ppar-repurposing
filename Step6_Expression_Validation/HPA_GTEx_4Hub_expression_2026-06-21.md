# HPA/GTEx 组织表达验证 — 4 Hub Genes

**日期**: 2026-06-21  
**来源**: Human Protein Atlas search API (`?format=json`), 整合 HPA + GTEx + FANTOM5  
**方法**: 文献校验标准 — 排除仅在无关组织（睾丸、脑等）高表达的基因，确认皮肤/脂肪中有表达

---

## 验证结果

| Gene | logFC (银屑病) | RNA 组织特异性 | 皮肤表达 | 脂肪表达 | 无关组织高表达 | 结论 |
|------|--------------|--------------|---------|---------|-------------|------|
| **ADIPOQ** | −1.449 | Group Enriched (Score=11) | ✅ 皮肤脂肪细胞 (nCPM 482.7) | ✅ 脂肪组织 (nTPM 843.2) | ❌ 未检测到睾丸/子-宫/脑 | ✅ **PASS** |
| **APOE** | −1.268 | Tissue Enhanced | ✅ 皮肤外根鞘细胞 | ✅ 肝脏 6533.9 nTPM | ⚠️ 脑 2714.9 nTPM（广谱表达） | ⚠️ **PASS (附注)** |
| **PLIN1** | −1.168 | Group Enriched (Score=7) | ✅ 皮肤脂肪细胞 (nCPM 879.6) | ✅ 脂肪组织 (nTPM 898.1) | ✅ 睾丸仅Sertoli细胞（体细胞）| ✅ **PASS** |
| **BCL2** | −1.100 | Low Specificity (广谱) | ✅ 广谱（Detected in all） | ✅ 脂肪细胞 (nCPM 1139.2) | ❌ 无睾丸/子宫富集 | ✅ **PASS** |

---

## 逐基因详细

### ADIPOQ (Adiponectin)
- **蛋白证据**: Evidence at protein level ✅
- **蛋白类别**: 分泌蛋白, 候选心血管病, 潜在药物靶点
- **皮肤表达**: HPA scRNA → Skin - Adipocytes (Skin) 显著富集
- **单细胞**: Adipocytes 482.7 nCPM（Cell type enriched）
- **排除组织**: 睾丸未检出, 脑未检出, 免疫细胞未检出
- **判断**: ✅ 皮肤皮下脂肪中明确表达。无无关组织污染。**通过**。

### APOE (Apolipoprotein E)
- **蛋白证据**: Evidence at protein level ✅
- **蛋白类别**: 脂代谢, 胆固醇代谢, 分泌蛋白
- **皮肤表达**: HPA scRNA → Skin - Outer root sheath cells 富集 ✅
- **问题**: 脑 (2714.9 nTPM)、肾上腺 (3391.3) 高表达。广谱表达基因，非皮肤/脂肪专一。
- **单细胞**: Hepatocytes 9962.2, Melanocytes 3427.4, Macrophages 1733.0
- **判断**: ⚠️ 皮肤有表达，但脑/肝表达远高于皮肤。解释为"系统脂蛋白代谢基因，其银屑病效应通过循环脂质→皮肤炎症轴介导，非局部靶组织必需专一表达"。

### PLIN1 (Perilipin 1)
- **蛋白证据**: Evidence at protein level ✅
- **蛋白类别**: 脂代谢, 质膜蛋白
- **皮肤表达**: HPA scRNA → Skin - Adipocytes (Skin) ✅
- **单细胞**: Adipocytes 879.6 nCPM
- **RNA**: 脂肪 898.1 nTPM, 睾丸 12.6 nTPM（低表达）
- **蛋白 IHC — 睾丸细胞类型分解**（HPA HTML tooltip 直接解析）:

| 睾丸细胞类型 | IHC 染色 | 细胞属性 |
|------------|---------|---------|
| Spermatogonia cells | Not detected | 生殖干细胞 |
| Preleptotene spermatocytes | Not detected | 生殖细胞 |
| Pachytene spermatocytes | Not detected | 生殖细胞 |
| Round/early spermatids | Not detected | 生殖细胞 |
| Elongated/late spermatids | Not detected | 生殖细胞 |
| Leydig cells | Not detected | 间质内分泌细胞 |
| Peritubular cells | Not detected | 管周肌样细胞 |
| **Sertoli cells** | **Medium** | **体细胞（支持/滋养）** |

- **判断**: ✅ 所有生殖细胞系 Not detected。Sertoli 细胞为体细胞支持细胞，不传递遗传物质。PLIN1 靶向药物生殖毒性风险极低。**通过**。

### BCL2 (BCL2 Apoptosis Regulator)
- **蛋白证据**: Evidence at protein level, FDA approved drug targets ✅
- **蛋白类别**: 凋亡/自噬, FDA 批准靶点, 原癌基因
- **皮肤表达**: 广谱 (Detected in all) ✅
- **单细胞**: Adipocytes 1139.2, B-cells 739.7, T-cells 738.2, Melanocytes 1004.4
- **排除组织**: 睾丸/子宫/脑无富集
- **判断**: ✅ 皮肤广谱表达。免疫细胞（B/T）表达支持银屑病免疫机制关联。

---

## 汇总

| | ADIPOQ | APOE | PLIN1 | BCL2 |
|---|---|---|---|---|
| 皮肤 RNA | ✅ Adipocytes | ✅ Outer root sheath | ✅ Adipocytes | ✅ 广谱 |
| 脂肪 RNA | ✅ 843.2 nTPM | ✅ | ✅ 898.1 nTPM | ✅ 1139.2 nCPM |
| 蛋白证据 | ✅ | ✅ | ✅ | ✅ FDA靶点 |
| 无关组织排除 | ✅ 通过 | ⚠️ 脑 (2714.9) | ⚠️ 睾丸精细胞 | ✅ 通过 |
| 最终判断 | **PASS** | **PASS (附说明)** | **PASS** ✅ | **PASS** |

### 已完成
- **PLIN1** ✅: HPA 蛋白 IHC 完整解析 — 睾丸仅 Sertoli 细胞（体细胞）Medium 染色，所有生殖细胞 Not detected。**通过**。原始数据归档：Human Protein Atlas 官网 PLIN1 页面 (https://www.proteinatlas.org/)
- **APOE**: 审稿人可能质疑脑高表达。需要在讨论中说明"银屑病中 APOE 效应通过循环脂蛋白→皮肤巨噬细胞炎症轴而非中枢神经系统介导"。

### 不需要排除的基因
4 个 Hub 基因全部有皮肤/脂肪表达证据，无需要硬排除的（如 CD36 旧管线中"睾丸 1439.4 nTPM vs 脂肪 688.2"那种压倒性差异）。PLIN1 是唯一有睾丸 co-expression 的，但睾丸也在皮肤脂肪细胞富集 → 不是"仅睾丸"的排除场景。
