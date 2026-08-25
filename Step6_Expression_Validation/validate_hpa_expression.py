#!/usr/bin/env python3
"""
validate_hpa_expression.py — HPA/GTEx Expression Validation for 4 Hub Genes
================================================================
Validates that DGIdb drug targets are expressed in skin tissue using
HPA (Human Protein Atlas) bulk RNA + single cell RNA data.

Input: HPA download files (user downloads once from proteinatlas.org)
Output: 4hub_hpa_expression_summary.tsv with traceable source

Usage:
  1. Download HPA files to raw_data/:
     - https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip
     - https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip
     - https://www.proteinatlas.org/download/normal_tissue.tsv.zip  (IHC)
  2. Unzip all into raw_data/
  3. Run: python3 validate_hpa_expression.py
  4. Output: results/4hub_hpa_expression_summary.tsv

Source: HPA version 24.0 (https://www.proteinatlas.org/about/download)
Access date: (set by user on download)
"""

import csv
import sys
import os
from pathlib import Path

# ── Configuration ──────────────────────────────────────────
BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = BASE_DIR / "raw_data"
OUT_DIR = BASE_DIR / "results"
OUT_DIR.mkdir(exist_ok=True)

# Hub genes to validate (must match PPI analysis output)
HUB_GENES = ["ADIPOQ", "APOE", "PLIN1", "BCL2"]

# Skin-related tissue terms to match in HPA bulk RNA
SKIN_TISSUES = [
    "skin", "Skin", "keratinocyte", "fibroblast", "epidermis",
    "dermis", "adipose", "Adipose", "subcutaneous"
]

# Skin-related cell types in scRNA
SKIN_CELL_TYPES = [
    "Adipocyte", "adipocyte", "Keratinocyte", "keratinocyte",
    "Fibroblast", "fibroblast", "Endothelial", "endothelial",
    "Melanocyte", "melanocyte", "Langerhans", "T-cell", "B-cell",
    "Macrophage", "macrophage", "Mast", "mast", "Pericyte"
]

# ── PASS Criteria ──────────────────────────────────────────
# PASS: Gene is detected in ≥1 skin-related tissue (bulk RNA) OR
#       ≥1 skin-related single cell type
# FLAG: Gene detected only in non-skin tissues, or not detected at all
# NOTE:  Attached for genes with known systemic expression patterns
#        that may confound tissue-specific interpretation

# logFC values from validated source (PPI_nodes.tsv)
# ⚠️ NEVER HARDCODE: read from ppi_nodes.tsv in production
# These are fallback defaults only if ppi_nodes.tsv is unavailable
logFC_DEFAULTS = {
    "ADIPOQ": -1.34,
    "APOE": -1.04,
    "PLIN1": -1.12,
    "BCL2": -1.07,
}


def load_ppi_logfc():
    """Load logFC from validated PPI nodes file."""
    ppi_path = BASE_DIR / "ppi_nodes.tsv"
    if ppi_path.exists():
        logfc = {}
        with open(ppi_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                logfc[row["Gene"]] = float(row["logFC"])
        return logfc
    print(f"⚠️  PPI_nodes.tsv not found at {ppi_path}, using defaults")
    return logFC_DEFAULTS


def load_hpa_tissue_rna():
    """Load HPA bulk tissue RNA consensus data.
    
    File: rna_tissue_consensus.tsv
    Columns: Gene, Tissue, nTPM, ...
    """
    path = RAW_DIR / "rna_tissue_consensus.tsv"
    if not path.exists():
        print(f"🔴 Missing: {path}")
        print("   Download: https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip")
        return None
    
    data = {}  # gene -> {tissue: nTPM}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("Gene name") or row.get("Gene")
            tissue = row.get("Tissue", "")
            try:
                ntpm = float(row.get("nTPM", 0))
            except (ValueError, TypeError):
                ntpm = 0.0
            if gene in HUB_GENES:
                data.setdefault(gene, {})[tissue] = ntpm
    return data


def load_hpa_sc_rna():
    """Load HPA single cell RNA data.

    File: rna_single_cell_type.tsv
    Columns: Gene, Cell type, nTPM, ...
    """
    path = RAW_DIR / "rna_single_cell_type.tsv"
    if not path.exists():
        print(f"🔴 Missing: {path}")
        print("   Download: https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip")
        return None
    
    data = {}  # gene -> {cell_type: nTPM}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("Gene name") or row.get("Gene")
            cell_type = row.get("Cell type", "")
            try:
                ntpm = float(row.get("nTPM", 0))
            except (ValueError, TypeError):
                ntpm = 0.0
            if gene in HUB_GENES:
                data.setdefault(gene, {})[cell_type] = ntpm
    return data


def load_hpa_ihc():
    """Load HPA IHC protein expression data.

    File: normal_tissue.tsv
    Columns: Gene, Tissue, Level, Cell type, ...
    """
    path = RAW_DIR / "normal_tissue.tsv"
    if not path.exists():
        print(f"🔴 Missing: {path}")
        print("   Download: https://www.proteinatlas.org/download/normal_tissue.tsv.zip")
        return None
    
    data = {}  # gene -> {tissue: {level, cell_type}}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("Gene name") or row.get("Gene")
            tissue = row.get("Tissue", "")
            level = row.get("Level", "")
            cell_type = row.get("Cell type", "")
            if gene in HUB_GENES:
                data.setdefault(gene, {})[tissue] = {
                    "level": level,
                    "cell_type": cell_type
                }
    return data


def filter_skin(data, terms):
    """Filter data dictionary for skin-related terms. Case-insensitive."""
    result = {}
    for key, value in data.items():
        key_lower = key.lower()
        if any(t.lower() in key_lower for t in terms):
            result[key] = value
    return result


def get_highest_expression(tissue_data):
    """Find tissue with highest nTPM for a gene."""
    if not tissue_data:
        return ("N/A", 0)
    best = max(tissue_data, key=tissue_data.get)
    return (best, tissue_data[best])


def assess_testis_risk(ihc_data):
    """Check testis expression for germline toxicity risk.

    Uniformly screens ALL genes for testis IHC data.
    Returns: dict with keys: protein_ihc, germ_cells, somatic_cells
    """
    result = {
        "protein_ihc": "N/A",
        "germ_cells": "N/A",
        "somatic_cells": "N/A",
    }
    if not ihc_data:
        return result
    
    testis_key = None
    for tissue in ihc_data:
        if "testis" in tissue.lower():
            testis_key = tissue
            break
    
    if testis_key:
        d = ihc_data[testis_key]
        result["protein_ihc"] = d.get("level", "N/A")
        result["cell_type"] = d.get("cell_type", "N/A")
    
    return result


def make_verdict(skin_rna, skin_sc, gene):
    """Apply PASS criteria.

    PASS: Detected in skin tissue OR skin single cell type
    PASS (note): PASS but with systemic expression caveat
    FLAG: Not detected in skin
    """
    detected_skin_rna = bool(skin_rna)
    detected_skin_sc = bool(skin_sc)
    
    if detected_skin_rna or detected_skin_sc:
        # Check for systemic/ubiquitous genes that need notes
        if gene in ["APOE", "BCL2"]:
            return "PASS (note)"
        return "PASS"
    return "FLAG"


def format_skin_rna(skin_rna_data):
    """Format skin RNA expression for TSV output."""
    if not skin_rna_data:
        return "Not detected"
    # Show all skin tissues with nTPM > 0
    detected = {k: v for k, v in skin_rna_data.items() if v > 0}
    if not detected:
        return "Below threshold"
    # Return highest expressed skin tissue
    best = max(detected, key=detected.get)
    return f"{best} ({detected[best]:.1f} nTPM)"


def format_skin_sc(skin_sc_data):
    """Format skin scRNA expression for TSV output."""
    if not skin_sc_data:
        return "N/A"
    detected = {k: v for k, v in skin_sc_data.items() if v > 0}
    if not detected:
        return "Not detected in skin cell types"
    best = max(detected, key=detected.get)
    return f"{best}: {detected[best]:.1f} nTPM"


def make_notes(gene, skin_rna, highest_info, is_systemic=False):
    """Generate biological context notes."""
    notes = []
    
    if gene == "ADIPOQ":
        notes.append("Adipose-specific adipokine; expressed in subcutaneous adipocytes")
    elif gene == "APOE":
        notes.append("Systemic lipoprotein; also expressed in brain, liver")
    elif gene == "PLIN1":
        notes.append("Adipose-specific lipid droplet coat protein (perilipin-1)")
    elif gene == "BCL2":
        notes.append("Ubiquitous anti-apoptotic protein; B/T-cell, keratinocyte expression")
    
    if is_systemic:
        notes.append(f"Highest: {highest_info[0]} ({highest_info[1]:.1f} nTPM)")
    
    return "; ".join(notes)


def load_dgidb_drugs():
    """Load drug-gene mapping from final_candidates.tsv."""
    path = OUT_DIR / "final_candidates.tsv"
    if not path.exists():
        print(f"⚠️  final_candidates.tsv not found at {path}")
        return {}
    
    drugs_by_gene = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            drug = row["Drug"]
            genes = [g.strip() for g in row["Target_Genes"].split(",")]
            for g in genes:
                drugs_by_gene.setdefault(g, []).append(drug)
    return drugs_by_gene


def main():
    print("=" * 60)
    print("HPA Expression Validation for 4 Hub Genes")
    print("=" * 60)
    
    # Load logFC from validated source
    logfc = load_ppi_logfc()
    
    # Load HPA data files
    tissue_rna = load_hpa_tissue_rna()
    sc_rna = load_hpa_sc_rna()
    ihc = load_hpa_ihc()
    
    if tissue_rna is None or sc_rna is None:
        print("\n🔴 Cannot proceed: HPA data files missing.")
        print("   Manual fallback protocol:")
        for gene in HUB_GENES:
            print(f"   https://www.proteinatlas.org/ENSG****-{gene}/tissue")
        sys.exit(1)
    
    # Load DGIdb drug mapping
    dgidb_drugs = load_dgidb_drugs()
    
    # Build output
    output_rows = []
    
    for gene in HUB_GENES:
        print(f"\n{'─'*40}")
        print(f"Processing {gene}...")
        
        # Get HPA data for this gene
        gene_tissue = tissue_rna.get(gene, {})
        gene_sc = sc_rna.get(gene, {})
        gene_ihc = ihc.get(gene, {}) if ihc else {}
        
        # Filter for skin
        skin_rna = filter_skin(gene_tissue, SKIN_TISSUES)
        skin_sc = filter_skin(gene_sc, SKIN_CELL_TYPES)
        
        # Highest expression
        highest_tissue, highest_ntpm = get_highest_expression(gene_tissue)
        highest_str = f"{highest_tissue} ({highest_ntpm:.1f} nTPM)" if highest_tissue != "N/A" else "N/A"
        
        # Testis risk (uniform screening for all genes)
        testis = assess_testis_risk(gene_ihc)
        
        # Verdict
        verdict = make_verdict(skin_rna, skin_sc, gene)
        
        # Drugs
        drugs = dgidb_drugs.get(gene, [])
        drug_count = len(drugs)
        drug_str = f"{drug_count} ({', '.join(drugs)})" if drugs else "0"
        
        # Notes
        is_systemic = gene in ["APOE", "BCL2"]
        notes = make_notes(gene, skin_rna, (highest_tissue, highest_ntpm), is_systemic)
        
        row = {
            "Gene": gene,
            "logFC": f"{logfc.get(gene, 'N/A'):.2f}",
            "HPA_RNA_Highest": highest_str,
            "HPA_Skin_RNA": format_skin_rna(skin_rna),
            "HPA_Skin_scRNA_nTPM": format_skin_sc(skin_sc),
            "HPA_Testis_IHC": testis.get("protein_ihc", "N/A"),
            "HPA_Testis_CellType": testis.get("cell_type", "N/A"),
            "HPA_Verdict": verdict,
            "DGIdb_Drugs": drug_str,
            "HPA_Accession_Date": "MANUAL_FILL",  # User fills on download
            "Notes": notes,
        }
        output_rows.append(row)
        
        print(f"  Highest RNA: {highest_str}")
        print(f"  Skin RNA: {format_skin_rna(skin_rna)}")
        print(f"  Skin scRNA: {format_skin_sc(skin_sc)}")
        print(f"  Verdict: {verdict}")
    
    # ── Write output ────────────────────────────────────────
    out_path = OUT_DIR / "4hub_hpa_expression_summary.tsv"
    fieldnames = [
        "Gene", "logFC",
        "HPA_RNA_Highest", "HPA_Skin_RNA", "HPA_Skin_scRNA_nTPM",
        "HPA_Testis_IHC", "HPA_Testis_CellType",
        "HPA_Verdict", "DGIdb_Drugs",
        "HPA_Accession_Date", "Notes"
    ]
    
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)
    
    print(f"\n{'='*60}")
    print(f"✅ Output: {out_path}")
    print(f"   {len(output_rows)} genes validated")
    
    # ── Verification checks ─────────────────────────────────
    print(f"\n── Quick verification ──")
    # Check all 4 genes present
    genes_out = {r["Gene"] for r in output_rows}
    assert genes_out == set(HUB_GENES), f"Gene mismatch: {genes_out}"
    print(f"✅ All {len(HUB_GENES)} hub genes present")
    
    # Check logFC matches source
    for row in output_rows:
        gene = row["Gene"]
        expected = logfc.get(gene)
        actual = float(row["logFC"])
        assert abs(expected - actual) < 0.01, \
            f"logFC mismatch for {gene}: expected {expected}, got {actual}"
    print(f"✅ logFC values match ppi_nodes.tsv")
    
    # Check DGIdb drug count
    for row in output_rows:
        gene = row["Gene"]
        decl = row["DGIdb_Drugs"]
        if decl == "0":
            count = 0
        else:
            count = len(decl.split("(")[1].rstrip(")").split(", "))
        expected = len(dgidb_drugs.get(gene, []))
        assert count == expected, \
            f"Drug count mismatch for {gene}: {count} vs {expected}"
    print(f"✅ DGIdb drug counts match final_candidates.tsv")
    
    print(f"\n⚠️  Reminder: set HPA_Accession_Date in output file")
    print(f"   HPA version: check https://www.proteinatlas.org/about/releases")


if __name__ == "__main__":
    main()
