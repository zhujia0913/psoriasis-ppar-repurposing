#!/usr/bin/env python3
"""
Batch docking of all 10 DGIdb candidate drugs against PPARα (PDB: 3VI8)
using AutoDock Vina v1.2.5 via the docking conda environment.
"""

import subprocess
import sys
import os
import shutil
from pathlib import Path

BASE = Path(__file__).resolve().parent
RECEPTOR = BASE / "receptors/PPARa_receptor.pdbqt"
RESULTS = BASE / "results"
LIGANDS = BASE / "ligands"
PYTHON = "python3"   # expects the 'docking' conda env (vina, rdkit, meeko) active
VINA = "vina"        # AutoDock Vina v1.2.x on PATH

# All 10 candidate drugs with PubChem CIDs
DRUGS = [
    ("spironolactone",    5833),
    ("fenofibrate",       3339),
    ("fluvastatin",       446155),
    ("irbesartan",        "irbesartan"),  # will try SDF download
    ("pravastatin",       54687),
    ("edaravone",         4021),
    ("isotretinoin",      5282379),
    ("ursodiol",          31401),
    ("erythromycin",      12560),
    ("pentoxifylline",    4740),
]

# PubChem CID lookup dict
PUBCHEM_CIDS = {
    "spironolactone": 5833,
    "fenofibrate": 3339,
    "fluvastatin": 446155,
    "irbesartan": 3742,
    "pravastatin": 54687,
    "edaravone": 4021,
    "isotretinoin": 5282379,
    "ursodiol": 31401,
    "erythromycin": 12560,
    "pentoxifylline": 4740,
}

def download_sdf(name, cid, dest):
    """Download SDF from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/SDF?record_type=3d"
    subprocess.run(["curl", "-s", "-L", url, "-o", str(dest)], check=True, timeout=30)
    size = dest.stat().st_size
    if size < 100:
        raise RuntimeError(f"SDF download too small ({size} bytes)")
    print(f"  Downloaded {name} SDF: {size} bytes")

def sdf_to_pdbqt(name, sdf_path, pdbqt_path):
    """Convert SDF to PDBQT using the docking env's meeko/rdkit."""
    script = f"""
import sys
    # rdkit/meeko are expected to be importable in the active (docking) environment
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

mol = None
supplier = Chem.SDMolSupplier('{sdf_path}')
for m in supplier:
    if m is not None:
        mol = m
        break

if mol is None:
    raise ValueError("Failed to load molecule")

mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)

preparator = MoleculePreparation()
preparator.prepare(mol)
preparator.write_pdbqt_file('{pdbqt_path}')
print(f"Written: {pdbqt_path}")
"""
    subprocess.run([str(PYTHON), "-c", script], check=True, timeout=60)

def run_docking(name):
    """Run AutoDock Vina docking for a single drug."""
    ligand_pdbqt = LIGANDS / f"{name}.pdbqt"
    out_pdbqt = RESULTS / f"{name}_docked.pdbqt"
    log_file = RESULTS / f"{name}_docking.log"
    
    # Write vina config
    config = f"""receptor = {RECEPTOR}
ligand = {ligand_pdbqt}
center_x = 11.060
center_y = 4.634
center_z = -7.613
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 32
num_modes = 10
out = {out_pdbqt}
"""
    config_path = BASE / "scripts" / f"vina_config_{name}.txt"
    config_path.write_text(config)
    
    # Run docking
    result = subprocess.run(
        [str(VINA), "--config", str(config_path)],
        capture_output=True, text=True, timeout=120
    )
    
    # Save log
    log_file.write_text(result.stdout + "\n" + result.stderr)
    
    # Parse best affinity
    for line in result.stdout.split("\n"):
        if line.strip().startswith("1 ") and "kcal" in line:
            parts = line.split()
            if len(parts) >= 2:
                return float(parts[1])
    
    return None

def main():
    skipped = []
    
    for name, cid_str in DRUGS:
        cid = PUBCHEM_CIDS.get(name)
        if cid is None:
            print(f"SKIP {name}: no CID")
            continue
        
        sdf_path = LIGANDS / f"{name}.sdf"
        pdbqt_path = LIGANDS / f"{name}.pdbqt"
        
        # Download SDF if not exists
        if not sdf_path.exists():
            try:
                download_sdf(name, cid, sdf_path)
            except Exception as e:
                print(f"FAIL {name} download: {e}")
                skipped.append(name)
                continue
        
        # Convert to PDBQT if not exists
        if not pdbqt_path.exists():
            try:
                sdf_to_pdbqt(name, sdf_path, pdbqt_path)
            except Exception as e:
                print(f"FAIL {name} prep: {e}")
                skipped.append(name)
                continue
        
        # Run docking if not already done
        log_file = RESULTS / f"{name}_docking.log"
        if log_file.exists():
            print(f"SKIP {name}: already docked")
            continue
        
        try:
            affinity = run_docking(name)
            if affinity is not None:
                print(f"DONE {name}: affinity = {affinity:.2f} kcal/mol")
        except Exception as e:
            print(f"FAIL {name} docking: {e}")
            skipped.append(name)
    
    print(f"\n=== SUMMARY ===")
    print(f"Skipped: {skipped}")
    
    # Parse all results
    print("\n=== ALL DOCKING RESULTS ===")
    for name, _ in DRUGS:
        log_file = RESULTS / f"{name}_docking.log"
        if not log_file.exists():
            print(f"{name}: NO LOG")
            continue
        text = log_file.read_text()
        for line in text.split("\n"):
            if line.strip().startswith("1 ") and "kcal" in line:
                parts = line.split()
                affinity = parts[1]
                print(f"{name}: {affinity} kcal/mol")
                break
        else:
            print(f"{name}: NO RESULT")

if __name__ == "__main__":
    main()
