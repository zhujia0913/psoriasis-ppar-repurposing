#!/usr/bin/env python3
"""
Convert SDF to PDBQT for AutoDock Vina using RDKit.
Handles partial charges and atom types for docking.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

def sdf_to_pdbqt(sdf_path, pdbqt_path):
    """Convert SDF to PDBQT format."""
    supplier = Chem.SDMolSupplier(sdf_path)
    mol = next(supplier)
    
    if mol is None:
        raise ValueError(f"Failed to load molecule from {sdf_path}")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates if not present
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
    
    # Write to PDB first
    pdb_path = pdbqt_path.replace('.pdbqt', '.pdb')
    Chem.MolToPDBFile(mol, pdb_path)
    
    # Convert PDB to PDBQT using Open Babel if available, otherwise use simple conversion
    try:
        import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb", "pdbqt")
        obmol = openbabel.OBMol()
        obConversion.ReadFile(obmol, pdb_path)
        obConversion.WriteFile(obmol, pdbqt_path)
        print(f"Converted {sdf_path} -> {pdbqt_path} (via Open Babel)")
    except ImportError:
        # Simple PDB to PDBQT conversion for docking
        # This is a simplified version - proper conversion needs partial charges
        pdb_to_pdbqt_simple(pdb_path, pdbqt_path)
        print(f"Converted {sdf_path} -> {pdbqt_path} (simple conversion)")
    
    # Clean up intermediate PDB
    if os.path.exists(pdb_path):
        os.remove(pdb_path)

def pdb_to_pdbqt_simple(pdb_path, pdbqt_path):
    """Simple PDB to PDBQT conversion without Open Babel."""
    with open(pdb_path, 'r') as f:
        pdb_lines = f.readlines()
    
    pdbqt_lines = []
    pdbqt_lines.append("REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.0\n")
    pdbqt_lines.append("REMARK  99 CREATED BY RDKIT SDF TO PDBQT CONVERTER\n")
    
    atom_num = 1
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Parse PDB ATOM line
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21:22].strip()
            res_num = line[22:26].strip()
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
            
            # Determine atom type (simplified)
            atom_type = atom_name[0] if atom_name else 'C'
            if atom_type in ['C', 'N', 'O', 'S', 'P', 'H', 'F', 'Cl', 'Br', 'I']:
                pass
            else:
                atom_type = 'A'  # Non-polar hydrogen or other
            
            # Create PDBQT ATOM line
            pdbqt_line = f"ATOM  {atom_num:5d}  {atom_name:4s} {res_name:3s} {chain:1s}{res_num:4s}    {x:8s}{y:8s}{z:8s}  1.00  0.00    {0.0:6.3f} {atom_type:2s}\n"
            pdbqt_lines.append(pdbqt_line)
            atom_num += 1
    
    pdbqt_lines.append("TER\n")
    pdbqt_lines.append("ENDMDL\n")
    
    with open(pdbqt_path, 'w') as f:
        f.writelines(pdbqt_lines)

if __name__ == "__main__":
    from pathlib import Path
    HERE = Path(__file__).resolve().parent
    ligands = HERE / "ligands"
    sdf_path = str(ligands / "spironolactone.sdf")
    pdbqt_path = str(ligands / "spironolactone.pdbqt")
    sdf_to_pdbqt(sdf_path, pdbqt_path)
    print(f"Output: {pdbqt_path}")
