# psoriasis-ppar-repurposing

Analysis code and pipelines accompanying the manuscript
*"Integrated transcriptomic and network analysis of peroxisome
proliferator-activated receptor target genes identifies repurposable drug
candidates in psoriasis"*.

> **Journal-independent repository.** This repository contains only the analysis
> code and the underlying numerical results; it is **not tied to any target
> journal**. Journal-specific packaging — section headings, supplementary-file
> names and numbering, reference style, figure formats, and the cover letter —
> is maintained separately in the submission package and is adapted if the
> manuscript is submitted elsewhere. The datasets, statistics, tables, and
> figure content generated here are identical across submissions: a change of
> target journal never requires re-running any analysis in this repository.

> **Correction notice (8 November 2026) — Step 9 / CIBERSORT.** The legacy
> `Step9_CIBERSORT/run_cibersort_analysis.R` contained two implementation bugs:
> (i) `as.vector(cor_hub)` flattened a 4 × 22 matrix column-wise while gene/cell
> labels were generated row-wise, systematically mismatching Rho and P values;
> (ii) `cor.test(method = "spearman")` exact p-values underflowed to exactly 0
> for strong correlations. Both are fixed in
> `Step9_CIBERSORT/CORRECTED_20261108/scripts/` (row-wise flattening via
> `as.vector(t(...))` and a t-approximation `2*pt(-abs(tstat), n-2)`). The
> corrected results (`CORRECTED_20261108/results/hub_immune_correlation_cibersort.tsv`,
> 41 significant pairs) are authoritative and match the manuscript, its
> hub–immune supplementary table, and Figs 12–13. The legacy script is retained for
> provenance only and must not be used for reproduction.

This repository provides the **complete, reproducible R/Python analysis scripts
for all ten analytical steps** that generate the findings and analytic tables
reported in the paper — from raw GEO retrieval through differential expression,
enrichment, network/hub-gene identification, immune deconvolution, drug
repositioning, and molecular docking. All scripts resolve paths **relative to
their own location**, so the repository is self-contained: clone it, place the
required public input files (see below) in the indicated folders, and run.

> **Figure numbering (manuscript):** Figure 14 = molecular docking of
> fenofibrate to PPARα; Figure 15 = molecular docking of pravastatin to PPARα.
> The earlier DGIdb drug–gene *network* figure is **not** included as a numbered
> figure; the DGIdb interaction results are reported in Table 2 and the
> accompanying drug-screening supplementary table.

---

## Repository structure

```
psoriasis-ppar-repurposing/
├── README.md
├── .gitignore
├── Step1_DEG/                  GEO retrieval, RMA, ComBat, limma DEG calling (Fig 2–4, 759 DEGs)
│   ├── deg_standard_pipeline_v2_AnnoProbe.R   GEO → RMA → ComBat → limma (→ expression_combat_corrected.csv, DEG lists)
│   └── pca_combat_v2_systematic.R_DEPRECATED   ⚠ OLD: synthetic pre-ComBat coords — DO NOT USE
├── Step2_PPAR_Overlap/          PPARgene ∩ DEG → 18-gene set (Fig 5)
│   └── extract_ppargene.py                    intersect PPARgene targets with DEGs → 18 genes
├── Step3_Enrichment/           GO & KEGG enrichment of 18 genes (Fig 6, 7)
│   └── run_enrichment_18_final.R              core 18-gene enrichment → data/*.tsv
├── Step4_PPI_MCODE/            STRING PPI + MCODE module (Fig 8, 9)
│   └── mcode_cytohubba.py                      PPI construction + MCODE module detection
├── Step5_Hub_CytoHubba/        CytoHubba rank-sum hub genes (Fig 10)
│   └── plot_ranksum_hub.py                      Fig 10 hub-gene rank-sum plot
├── Step6_Expression_Validation/  HPA/GTEx tissue expression (Table 1)
│   ├── run_tau_dgidb_integration.R             TAU specificity index + DGIdb integration
│   └── validate_hpa_expression.py             HPA validation
├── Step7_DGIdb_DrugRepurposing/  DGIdb querying + 5-layer filtering (Table 2 / drug-screening supplementary table)
│   └── DGIdb_v5_final_candidates.tsv           final 10 candidates
├── Step8_TRRUST_TF_Network/     TRRUST TF regulatory network (Fig 11)
│   ├── scripts/build_trrust_network.py         build TF–hub network from TRRUST
│   ├── data/                                   TRRUST edges/nodes (processed)
│   ├── raw_data/trrust_rawdata_human.tsv       TRRUST v2 raw dump
│   └── figures/                                network tables
├── Step9_CIBERSORT/            Immune-cell deconvolution (Fig 12, 13)
│   ├── CORRECTED_20261108/     authoritative Step 9 scripts + results (USE THIS)
│   └── run_cibersort_analysis.R               DEPRECATED buggy version - do not use
├── Step10_Docking/             Molecular docking (Fig 14 / 15)
│   ├── batch_dock.py                           AutoDock Vina docking of 10 candidates → PPARα LBD
│   └── sdf_to_pdbqt.py                         ligand prep (SDF → PDBQT)
└── Step11_External_Validation/ Independent-cohort validation (GSE30999; Fig 16)
    ├── scripts/                                00 fetch → 01 preprocess → 02 validate
    ├── data/                                   runtime GEO inputs (not committed)
    └── results/                                external-validation table + Fig 16
```

Each step reads inputs from and writes outputs to `data/` and/or `figures/`
sub-folders **next to the script**. Raster images are not committed (see
`.gitignore`); the small tabular intermediate products that support
reproducibility are committed in each step's `data/` folder.

**Note on figure-rendering scripts:** the 15 presentation-layer scripts that
draw Fig 2–15 (PCA / volcano / heatmap / Venn / GO / KEGG / PPI / violin /
Spearman-heatmap / PyMOL-render / composite scripts) are *not* included in this
repository. They are not required to reproduce the reported quantitative
results, which are produced by the analysis scripts listed above; they are
archived separately outside this repository for reference.

---

## Input data & how to obtain them

The scripts are self-contained *given the processed inputs below*. Raw public
datasets are **not** redistributed; obtain them from the cited resources and
place them in the indicated relative paths.

| Input | Source / how to obtain | Place at |
|-------|----------------------|----------|
| **Raw GEO CEL files** (GSE13355, GSE14905, GSE78097) | Download from GEO; Affymetrix HG-U133 Plus 2.0 / GPL570. 209 skin samples. | `Step1_DEG/data/raw/` (set `RAW_DIR` in `deg_standard_pipeline_v2_AnnoProbe.R`) |
| **ComBat-corrected expression matrix** (`expression_combat_corrected.csv`) | Generated by Step 1 (RMA → ComBat → limma). 68 MB; not committed. Regenerate via Step 1, or download from the Zenodo deposit cited in the manuscript (link and DOI to be added once the archive is created). | `Step1_DEG/expression_combat_corrected.csv` |
| **CIBERSORT LM22 signature matrix** (`LM22.txt`) | **License-restricted** — register at <https://cibersortx.stanford.edu/>. Cannot be redistributed in this repo. | `Step9_CIBERSORT/data/LM22.txt` (+ `CIBERSORT.R`) |
| **PPARα LBD structure** (PDB **3VI8**, 1.75 Å) | Download from RCSB PDB (<https://www.rcsb.org/structure/3VI8>). Prepare receptor (remove water/ligand, add polar hydrogens) → `PPARa_receptor.pdbqt`. | `Step10_Docking/receptors/` |
| **Ligand 3D structures** | Fetched automatically by `batch_dock.py` from PubChem (SDF, 3D) using the CIDs listed in the script (e.g. fenofibrate 3339, pravastatin 54687). | `Step10_Docking/ligands/` (created at runtime) |
| **DGIdb v5.0 drug–gene interactions** | Queried at <https://www.dgidb.org/>; raw response included as `dgidb_v5_4hub_raw_response.json`. | not required to re-run docking |
| **TRRUST v2 raw data** (`trrust_rawdata_human.tsv`) | Downloaded from TRRUST; included in `Step8_TRRUST_TF_Network/raw_data/`. | — |
| **GSE30999 processed files** (external validation cohort) | Download `E-GEOD-30999.processed.1.zip` (38 MB) + `E-GEOD-30999.sdrf.txt` from ArrayExpress / BioStudies (Affymetrix GPL570; 85 patients / 170 paired biopsies). Not redistributed. | `Step11_External_Validation/data/` (or set `GSE30999_DATA_DIR`; fetched automatically by `scripts/00_fetch_data.R`) |
| **Enrichment gene list (18 genes)** | **Embedded directly in `run_enrichment_18_final.R`**. | — |
| **GO / KEGG reference** (org.Hs.eg.db, KEGG) | Retrieved live by clusterProfiler / org.Hs.eg.db at runtime (internet required). | — |

**Databases used for interpretation (no local files needed):**
PPARgene (177 experimentally verified human PPAR target genes),
STRING v12 (PPI), HPA & GTEx (tissue expression),
TRRUST v2 (TF–target regulatory network), DGIdb v5.0 (drug–gene).

---

## How to reproduce

### Step 1 — DEG identification (Fig 2–4)
```bash
cd Step1_DEG
# place GEO CEL files in data/raw/ (see table above)
Rscript deg_standard_pipeline_v2_AnnoProbe.R   # → expression_combat_corrected.csv, deg_results_*.tsv
```

### Step 2 — PPARgene ∩ DEG → 18 genes (Fig 5)
```bash
cd Step2_PPAR_Overlap
python extract_ppargene.py                     # → PPARgene_DEG_overlap_18.tsv
```

### Step 3 — Enrichment (Fig 6, 7)
```bash
cd Step3_Enrichment
Rscript run_enrichment_18_final.R      # writes data/GO_*.tsv, KEGG_enrichment.tsv
```
Arial must be installed (macOS: `/System/Library/Fonts/Supplemental/Arial*.ttf`;
other OS: adjust the `font_add()` paths or install Arial).

### Step 4 — PPI + MCODE (Fig 8, 9)
```bash
cd Step4_PPI_MCODE
python mcode_cytohubba.py              # builds PPI, detects MCODE module
```

### Step 5 — Hub genes (Fig 10)
```bash
cd Step5_Hub_CytoHubba
python plot_ranksum_hub.py             # Fig 10
```

### Step 6 — Tissue expression validation (Table 1)
```bash
cd Step6_Expression_Validation
Rscript run_tau_dgidb_integration.R     # TAU specificity index + DGIdb integration → results/*.tsv
python validate_hpa_expression.py       # HPA expression summary → results/4hub_hpa_expression_summary.tsv
```

### Step 7 — DGIdb drug repositioning (Table 2 / drug-screening supplementary table)
```bash
cd Step7_DGIdb_DrugRepurposing
# DGIdb_v5_final_candidates.tsv = final 10 candidates (Table 2 / drug-screening supplementary table)
# Raw API responses and the 5-layer filtering pipeline are documented in
# METHODOLOGY.md; the DGIdb query script is not included in this repository.
```

### Step 8 — TRRUST TF network (Fig 11)
```bash
cd Step8_TRRUST_TF_Network
python scripts/build_trrust_network.py # Fig 11 (reads raw_data/, writes data/ & figures/)
```

### Step 9 — CIBERSORT (Fig 12, 13) — USE THE CORRECTED VERSION
```bash
cd Step9_CIBERSORT/CORRECTED_20261108
# place raw_data/expression_pre_combat.csv (68 MB; regenerate via Step 1),
# data/sample_info.tsv (included) and data/LM22.txt (license-restricted)
Rscript scripts/run_cibersort_analysis_rerun.R   # writes results/ + figures/
```
The legacy `Step9_CIBERSORT/run_cibersort_analysis.R` is **deprecated**; see the
Correction notice at the top of this README.

### Step 10 — Molecular docking (Fig 14, 15)
```bash
cd Step10_Docking
# activate the conda env that provides AutoDock Vina, RDKit, Meeko, PyMOL
conda activate docking

python batch_dock.py                      # docks 10 ligands → results/*_docked.pdbqt + logs
python sdf_to_pdbqt.py                    # (helper) convert a single SDF → PDBQT
```

### Step 11 — External validation in an independent cohort (GSE30999 / Fig 16)
```bash
cd Step11_External_Validation
Rscript scripts/00_fetch_data.R          # ArrayExpress E-GEOD-30999 (38 MB; internet)
Rscript scripts/01_preprocess.R          # MaxMean probe selection + 81 paired samples
Rscript scripts/02_external_validation.R # paired Wilcoxon + DeLong AUC → results/ + figures/
```
The pre-specified 18-gene and 4-hub inputs, plus the 18 MaxMean representative
probes, are in `scripts/input/`. `results/external_validation_table.xlsx` is a
byte-identical reference copy of the supplementary external-validation table.
All reported numbers are reproducible end to end from the code above; see
`Step11_External_Validation/README.md` and
`Step11_External_Validation/results/reproducibility_check.tsv` (R vs independent
Python recomputation, 18/18 rows identical).

---

## Requirements

- **R** (≥ 4.0): affy/oligo, sva, limma, AnnoProbe, ggplot2, dplyr, cowplot,
  showtext, grid, gtable, stringr, scales, patchwork, clusterProfiler,
  org.Hs.eg.db, enrichplot, pheatmap, reshape2, tidyr, ggpubr, corrplot,
  CIBERSORT (LM22), GEOquery, hgu133plus2.db, pROC. *(Steps 1, 3, 6, 9, 11.)*
- **Python** (≥ 3.8): matplotlib, pandas, numpy, networkx;
  PyMOL (rendering) and a docking environment with **AutoDock Vina ≥ 1.2**,
  **RDKit**, **Meeko** (e.g. a `conda` env named `docking`).
  *(Steps 2, 4, 5, 8, 10.)*
- Step 4 MCODE module detection uses the MCODE algorithm (igraph implementation
  or the Cytoscape MCODE plugin); the module intermediates are committed, so the
  MCODE result can be regenerated without Cytoscape.

## Output format

All raster figures are exported as **TIFF, 300–600 DPI, LZW compression,
flattened RGB, < 10 MB per file**. Combination figures (e.g. docking panels) use
300 DPI; line/dot plots use up to 600 DPI. Figure fonts use **Arial** (8–12 pt).
These are objective, journal-independent export settings that satisfy the figure
requirements of most biomedical journals; if a target journal specifies a
different format or resolution, re-export from the scripts rather than changing
any underlying data.

## Data Availability

All analysis code in this repository is released under the MIT license and is
**publicly available**. The key intermediate analytic tables that support the
manuscript's findings — the 759-gene DEG list, the 18-gene set, CIBERSORT immune
scores, DGIdb candidates, TRRUST edges, and the GSE30999 external-validation
statistics — are provided both in each step's
`data/` folder and as supplementary tables accompanying the manuscript. The
68 MB Pre-ComBat expression matrix required by Step 9 is regenerated by Step 1
or obtained from the Zenodo deposit cited in the manuscript (a version-stable
archive with a DOI is planned; the link will be added here once available).

## License

MIT (code). Supplementary analytic tables accompanying the manuscript are
distributed under the same terms.
