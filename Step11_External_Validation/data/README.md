# data/

Runtime inputs for the GSE30999 external-validation stage. **No GEO data are
committed to this repository**; the `.rds` and series-matrix files are listed
in `.gitignore`.

From `Step11_External_Validation/`, run:

```bash
Rscript scripts/00_fetch_data.R
```

This downloads (or reuses) `GSE30999_series_matrix.txt.gz` and writes:

- `gse30999_expression.rds` — probes × samples expression matrix (log2, RMA)
- `gse30999_pheno.rds` — sample metadata

`scripts/01_preprocess.R` then writes `gse30999_processed.rds` (gene × sample)
and `gse30999_pairs.rds` (81 paired lesional / non-lesional samples).

## Source

- GEO accession: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30999>
- Platform: Affymetrix HG-U133 Plus 2.0 (GPL570)
- 85 patients; 170 paired biopsies; 81 pairs after QC.
