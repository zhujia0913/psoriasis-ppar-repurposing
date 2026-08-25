# psoriasis-ppar-repurposing

Analysis code and pipelines accompanying the manuscript
*"A PPAR-centric framework for lipid–immune crosstalk in psoriasis and
repositioned therapeutics"* (submitted to PLOS ONE).

This repository provides the **complete, reproducible R/Python scripts for all
ten analytical steps** that generate the figures and analytic tables reported in
the paper — from raw GEO retrieval through differential expression, enrichment,
network/hub-gene identification, immune deconvolution, drug repositioning, and
molecular docking. All scripts resolve paths **relative to their own location**,
so the repository is self-contained: clone it, place the required public input
files (see below) in the indicated folders, and run.

> **Figure numbering (final manuscript):** Figure 14 = molecular docking of
> fenofibrate to PPARα; Figure 15 = molecular docking of pravastatin to PPARα.
> The earlier DGIdb drug–gene *network* figure is **not** included as a numbered
> figure; the DGIdb interaction results are reported in Table 2 / S2 Table.

---

## Repository structure

```
psoriasis-ppar-repurposing/
├── README.md
├── .gitignore
├── Step1_DEG/                  GEO retrieval, RMA, ComBat, limma DEG calling (Fig 2–4, 759 DEGs)
│   ├── deg_standard_pipeline_v2_AnnoProbe.R   GEO → RMA → ComBat → limma (→ expression_combat_corrected.csv, DEG lists)
│   ├── pca_combat_v2_systematic.R             Fig 2 PCA before/after ComBat
│   ├── volcano_deg.R                           Fig 3 volcano plot
│   └── heatmap_top100.R                        Fig 4 top-100 DEG heatmap
├── Step2_PPAR_Overlap/          PPARgene ∩ DEG → 18-gene set (Fig 5)
│   ├── extract_ppargene.py                    intersect PPARgene targets with DEGs → 18 genes
│   └── draw_venn.py                            Fig 5 Venn diagram
├── Step3_Enrichment/           GO & KEGG enrichment of 18 genes (Fig 6, 7)
│   ├── run_enrichment_18_final.R              core 18-gene enrichment → data/*.tsv
│   ├── run_3x2_reference_style.R              Fig 6 GO 3×2 panel (Arial, 3000 DPI)
│   ├── kegg_enrichment_plot.R                 Fig 7 KEGG plots
│   └── fix_dpi_tiff.py                         PNG → TIFF (LZW)
├── Step4_PPI_MCODE/            STRING PPI + MCODE module (Fig 8, 9)
│   ├── mcode_cytohubba.py                      PPI construction + MCODE module detection
│   ├── plot_ppi_network_hub.py                 Fig 8 PPI network
│   └── plot_ppi_module1.py                     Fig 9 MCODE module
├── Step5_Hub_CytoHubba/        CytoHubba rank-sum hub genes (Fig 10)
│   └── plot_ranksum_hub.py                      Fig 10 hub-gene rank-sum plot
├── Step6_Expression_Validation/  HPA/GTEx tissue expression (Table 1)
│   ├── run_tau_dgidb_integration.R             TAU specificity index + DGIdb integration
│   └── validate_hpa_expression.py             HPA validation
├── Step7_DGIdb_DrugRepurposing/  DGIdb querying + 5-layer filtering (Table 2 / S2)
│   └── DGIdb_v5_final_candidates.tsv           final 10 candidates
├── Step8_TRRUST_TF_Network/     TRRUST TF regulatory network (Fig 11)
│   ├── scripts/build_trrust_network.py         build TF–hub network from TRRUST
│   ├── data/                                   TRRUST edges/nodes (processed)
│   ├── raw_data/trrust_rawdata_human.tsv       TRRUST v2 raw dump
│   └── figures/                                network tables
├── Step9_CIBERSORT/            Immune-cell deconvolution (Fig 12, 13)
│   ├── run_cibersort_analysis.R               CIBERSORTx on psoriatic vs control
│   ├── plot_violin_literature_style.R          Fig 12 22-cell violin (Arial, 3000 DPI)
│   ├── redraw_heatmap_hub_immune.R             Fig 13 hub × immune Spearman heatmap
│   └── fix_dpi_tiff.py                         PNG → TIFF (LZW)
└── Step10_Docking/             Molecular docking (Fig 14 / 15)
    ├── batch_dock.py                           AutoDock Vina docking of 10 candidates → PPARα LBD
    ├── sdf_to_pdbqt.py                         ligand prep (SDF → PDBQT)
    ├── render_fenofibrate_3000dpi.py           PyMOL render of fenofibrate complex
    ├── render_pravastatin_3000dpi.py           PyMOL render of pravastatin complex
    └── composite_plos.py                        PLOS ONE–compliant TIFF composite (300 DPI, LZW)
```

Each step reads inputs from and writes outputs to `data/` and/or `figures/`
sub-folders **next to the script**. Raster images are not committed (see
`.gitignore`); the small tabular intermediate products that support
reproducibility are committed in each step's `data/` folder.

---

## Input data & how to obtain them

The scripts are self-contained *given the processed inputs below*. Raw public
datasets are **not** redistributed; obtain them from the cited resources and
place them in the indicated relative paths.

| Input | Source / how to obtain | Place at |
|-------|----------------------|----------|
| **Raw GEO CEL files** (GSE13355, GSE14905, GSE78097) | Download from GEO; Affymetrix HG-U133 Plus 2.0 / GPL570. 209 skin samples. | `Step1_DEG/data/raw/` (set `RAW_DIR` in `deg_standard_pipeline_v2_AnnoProbe.R`) |
| **ComBat-corrected expression matrix** (`expression_combat_corrected.csv`) | Generated by Step 1 (RMA → ComBat → limma). 68 MB; not committed. Regenerate via Step 1, or download from the Zenodo/Supporting Information deposit cited in the manuscript. | `Step1_DEG/expression_combat_corrected.csv` |
| **CIBERSORT LM22 signature matrix** (`LM22.txt`) | **License-restricted** — register at <https://cibersortx.stanford.edu/>. Cannot be redistributed in this repo. | `Step9_CIBERSORT/data/LM22.txt` (+ `CIBERSORT.R`) |
| **PPARα LBD structure** (PDB **3VI8**, 1.75 Å) | Download from RCSB PDB (<https://www.rcsb.org/structure/3VI8>). Prepare receptor (remove water/ligand, add polar hydrogens) → `PPARa_receptor.pdbqt`. | `Step10_Docking/receptors/` |
| **Ligand 3D structures** | Fetched automatically by `batch_dock.py` from PubChem (SDF, 3D) using the CIDs listed in the script (e.g. fenofibrate 3339, pravastatin 54687). | `Step10_Docking/ligands/` (created at runtime) |
| **DGIdb v5.0 drug–gene interactions** | Queried at <https://www.dgidb.org/>; raw response included as `dgidb_v5_4hub_raw_response.json`. | not required to re-run docking |
| **TRRUST v2 raw data** (`trrust_rawdata_human.tsv`) | Downloaded from TRRUST; included in `Step8_TRRUST_TF_Network/raw_data/`. | — |
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
Rscript pca_combat_v2_systematic.R             # Fig 2
Rscript volcano_deg.R                          # Fig 3
Rscript heatmap_top100.R                       # Fig 4
```

### Step 2 — PPARgene ∩ DEG → 18 genes (Fig 5)
```bash
cd Step2_PPAR_Overlap
python extract_ppargene.py                     # → PPARgene_DEG_overlap_18.tsv
python draw_venn.py                            # Fig 5
```

### Step 3 — Enrichment (Fig 6, 7)
```bash
cd Step3_Enrichment
Rscript run_enrichment_18_final.R      # writes data/GO_*.tsv, KEGG_enrichment.tsv
Rscript run_3x2_reference_style.R      # reads data/, writes figures/ (PNG+PDF+TIFF)
Rscript kegg_enrichment_plot.R         # reads data/, writes figures/
```
Arial must be installed (macOS: `/System/Library/Fonts/Supplemental/Arial*.ttf`;
other OS: adjust the `font_add()` paths or install Arial).

### Step 4 — PPI + MCODE (Fig 8, 9)
```bash
cd Step4_PPI_MCODE
python mcode_cytohubba.py              # builds PPI, detects MCODE module
python plot_ppi_network_hub.py         # Fig 8
python plot_ppi_module1.py             # Fig 9
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

### Step 7 — DGIdb drug repositioning (Table 2 / S2)
```bash
cd Step7_DGIdb_DrugRepurposing
# DGIdb_v5_final_candidates.tsv = final 10 candidates (Table 2 / S2)
# Raw API responses and the 5-layer filtering pipeline are documented in
# METHODOLOGY.md; the DGIdb query script is not included in this repository.
```

### Step 8 — TRRUST TF network (Fig 11)
```bash
cd Step8_TRRUST_TF_Network
python scripts/build_trrust_network.py # Fig 11 (reads raw_data/, writes data/ & figures/)
```

### Step 9 — CIBERSORT (Fig 12, 13)
```bash
cd Step9_CIBERSORT
# place data/expression_combat_corrected.csv and data/LM22.txt (see table above)
Rscript run_cibersort_analysis.R          # writes data/ + results/ + figures/
Rscript plot_violin_literature_style.R    # Fig 12 → figures/ (3000 DPI TIFF)
Rscript redraw_heatmap_hub_immune.R       # Fig 13 → figures/
```

### Step 10 — Molecular docking (Fig 14, 15)
```bash
cd Step10_Docking
# activate the conda env that provides AutoDock Vina, RDKit, Meeko, PyMOL
conda activate docking

python batch_dock.py                      # docks 10 ligands → results/*_docked.pdbqt + logs
python sdf_to_pdbqt.py                    # (helper) convert a single SDF → PDBQT
python render_fenofibrate_3000dpi.py      # PyMOL → figures/fenofibrate_panelA/B_3000dpi.png
python render_pravastatin_3000dpi.py      # PyMOL → figures/pravastatin_panelA/B_3000dpi.png
python composite_plos.py                  # figures/ → plos_tiff/Figure_14/15_*.tif (300 DPI, LZW)
```

---

## Requirements

- **R** (≥ 4.0): affy/oligo, sva, limma, AnnoProbe, ggplot2, dplyr, cowplot,
  showtext, grid, gtable, stringr, scales, patchwork, clusterProfiler,
  org.Hs.eg.db, enrichplot, pheatmap, reshape2, tidyr, ggpubr, corrplot,
  CIBERSORT (LM22). *(Steps 1, 3, 6, 9.)*
- **Python** (≥ 3.8): Pillow (TIFF/DPI), matplotlib, pandas, numpy, networkx;
  PyMOL (rendering) and a docking environment with **AutoDock Vina ≥ 1.2**,
  **RDKit**, **Meeko** (e.g. a `conda` env named `docking`).
  *(Steps 2, 4, 5, 7, 8, 10.)*
- Step 4 MCODE module detection uses the MCODE algorithm (igraph implementation
  or the Cytoscape MCODE plugin); the module intermediates are committed, so the
  plotting scripts run without Cytoscape.

## Output format

All raster figures meet the **PLOS ONE** specification: **TIFF, 300–600 DPI,
LZW compression, flattened RGB**. Combination figures (e.g. docking panels) use
300 DPI; line/dot plots use up to 600 DPI. Figure fonts use **Arial** (8–12 pt).

## Data Availability

All analysis code in this repository is released under the MIT license and is
**publicly available**. The key intermediate analytic tables that support the
manuscript's findings — the 759-gene DEG list, the 18-gene set, CIBERSORT immune
scores, DGIdb candidates, and TRRUST edges — are provided both in each step's
`data/` folder and as Supporting Information (S1–S3) accompanying the manuscript.
The 68 MB ComBat-corrected expression matrix is regenerated by Step 1 or obtained
from the Zenodo deposit cited in the manuscript. **Make the repository public (or
archive a release on Zenodo) before acceptance.**

## License

MIT (code). Supporting analytic tables (S1–S3) accompany the manuscript as
Supporting Information.
