# data/

Runtime inputs for the GSE30999 external-validation stage. **No source data
are committed to this repository**; the processed files and the `.rds`
objects are listed in `.gitignore`.

From `Step11_External_Validation/`, run:

```bash
Rscript scripts/00_fetch_data.R
```

This downloads (or reuses) the two ArrayExpress / BioStudies files for
`E-GEOD-30999` and verifies their sizes and the SDRF md5:

| File | Bytes | Notes |
|------|-------|-------|
| `E-GEOD-30999.processed.1.zip` | 39,798,039 | 170 `GSM*_sample_table.txt` files (`ID_REF / VALUE`), 54,675 probes, RMA log2 values |
| `E-GEOD-30999.sdrf.txt` | 109,410 | md5 `55d99d97c869f534a6f29ea6df85437e` |

Official source:
<https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/GEOD/E-GEOD-30999/>

If the files are already on disk, point the scripts at their directory:

```bash
GSE30999_DATA_DIR=/path/to/data Rscript scripts/00_fetch_data.R
GSE30999_DATA_DIR=/path/to/data Rscript scripts/01_preprocess.R
```

`scripts/01_preprocess.R` then writes `gse30999_processed.rds`
(gene × sample, MaxMean collapse) and `gse30999_pairs.rds` (81 paired
lesional / non-lesional samples). Both are git-ignored intermediate objects.

## Source

- GEO accession: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30999>
- ArrayExpress / BioStudies accession: `E-GEOD-30999`
- Platform: Affymetrix HG-U133 Plus 2.0 (GPL570)
- 85 patients; 170 paired biopsies; 81 pairs after QC.
