#!/usr/bin/env python3
"""
Pravastatin-PPARα 2-panel figure — Journal Quality (3000 DPI, Arial, Clean Aesthetics)
Panel A: Overview (cartoon + ligand + pocket label)
Panel B: Binding site close-up (sticks + dashes + residue labels)
Style: receptor=slateblue, ligand=firebrick, binding_site=lightsteelblue, Arial labels.
"""
import pymol
from pymol import cmd
import math
import os

pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

from pathlib import Path
HERE = Path(__file__).resolve().parent
RECEPTOR = str(HERE / "receptors" / "PPARa_clean.pdb")
LIGAND = str(HERE / "results" / "pravastatin_best.pdb")
OUT_DIR = str(HERE / "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ============ Panel A: Overview ============
cmd.reinitialize()
cmd.set_color("slateblue_custom", [0.30, 0.38, 0.55])
cmd.set_color("lightsteelblue_custom", [0.69, 0.77, 0.87])
cmd.load(RECEPTOR, "receptor")
cmd.load(LIGAND, "ligand")

# Premium white background, maximum quality
cmd.bg_color("white")
cmd.set("antialias", 4)
cmd.set("ray_trace_mode", 1)
cmd.set("orthoscopic", "on")
cmd.set("cartoon_smooth_loops", "off")
cmd.set("ray_shadows", "off")

# Receptor: elegant slate blue
cmd.hide("everything")
cmd.show("cartoon", "receptor")
cmd.color("slateblue_custom", "receptor")
cmd.set("cartoon_transparency", 0.25)
cmd.set("cartoon_fancy_helices", "on")
cmd.set("cartoon_highlight_color", "white")

# Ligand: rich firebrick
cmd.show("sticks", "ligand")
cmd.show("spheres", "ligand")
cmd.color("firebrick", "ligand")
cmd.set("stick_radius", 0.15)
cmd.set("sphere_scale", 0.18, "ligand")

# Calculate ligand center
model = cmd.get_model("ligand")
coords = [a.coord for a in model.atom]
cx = sum(c[0] for c in coords) / len(coords)
cy = sum(c[1] for c in coords) / len(coords)
cz = sum(c[2] for c in coords) / len(coords)

# Pocket label with Arial
label_pos = (cx + 20, cy + 20, cz + 8)
cmd.pseudoatom("pocket_label", pos=label_pos)
cmd.label("pocket_label", '"Binding pocket"')
try:
    cmd.set("label_font_name", "Arial", "pocket_label")
except:
    cmd.set("label_font_id", 7, "pocket_label")
cmd.set("label_size", 20, "pocket_label")
cmd.set("label_color", "black", "pocket_label")
cmd.set("label_outline_color", "white", "pocket_label")

# Elegant dashed arrow
def add_dashed_arrow(name, start, end, color=(0.8, 0.2, 0.2), dash_len=1.2, gap=0.6, radius=0.10):
    x1, y1, z1 = start
    x2, y2, z2 = end
    dx, dy, dz = x2-x1, y2-y1, z2-z1
    total = math.sqrt(dx*dx + dy*dy + dz*dz)
    if total == 0: return
    nx, ny, nz = dx/total, dy/total, dz/total
    obj = []
    t = 0.0
    while t < total:
        s0 = t
        s1 = min(t + dash_len, total)
        if s1 - s0 < 0.2: break
        obj.extend([
            25.0, color[0], color[1], color[2],
            x1+nx*s0, y1+ny*s0, z1+nz*s0,
            x1+nx*s1, y1+ny*s1, z1+nz*s1,
            radius,
            color[0], color[1], color[2],
            color[0], color[1], color[2],
            1.0, 1.0
        ])
        t += dash_len + gap
    cmd.load_cgo(obj, name)

add_dashed_arrow("arrow", label_pos, (cx, cy, cz), color=(0.8, 0.2, 0.2))

# Receptor label with Arial
label2_pos = (cx - 28, cy - 28, cz - 12)
cmd.pseudoatom("receptor_label", pos=label2_pos)
cmd.label("receptor_label", '"PPAR\u03b1 (PDB: 3VI8)"')
try:
    cmd.set("label_font_name", "Arial", "receptor_label")
except:
    cmd.set("label_font_id", 7, "receptor_label")
cmd.set("label_size", 20, "receptor_label")
cmd.set("label_color", "black", "receptor_label")
cmd.set("label_outline_color", "white", "receptor_label")

# Optimal view
cmd.center("ligand")
cmd.zoom("receptor", 14)

# 3000 DPI rendering
cmd.ray(3600, 2800)
panel_a = os.path.join(OUT_DIR, "pravastatin_panelA_3000dpi.png")
cmd.png(panel_a, dpi=3000, width=3600, height=2800)
print(f"Panel A done: {panel_a}")

# ============ Panel B: Binding site close-up ============
cmd.reinitialize()
cmd.set_color("slateblue_custom", [0.30, 0.38, 0.55])
cmd.set_color("lightsteelblue_custom", [0.69, 0.77, 0.87])
cmd.load(RECEPTOR, "receptor")
cmd.load(LIGAND, "ligand")

cmd.bg_color("white")
cmd.set("antialias", 4)
cmd.set("ray_trace_mode", 1)
cmd.set("orthoscopic", "on")
cmd.set("cartoon_smooth_loops", "off")
cmd.set("ray_shadows", "off")

# Binding site
cmd.hide("everything", "receptor")
cmd.select("binding_site", "byres (receptor within 4.0 of ligand)")
cmd.show("sticks", "binding_site")
cmd.color("lightsteelblue_custom", "binding_site")

# Ligand
cmd.show("sticks", "ligand")
cmd.show("spheres", "ligand")
cmd.color("firebrick", "ligand")
cmd.set("stick_radius", 0.15)
cmd.set_bond("stick_radius", 0.15, "ligand", "ligand")
cmd.set_bond("stick_radius", 0.12, "binding_site", "binding_site")
cmd.set("sphere_scale", 0.18, "ligand")

# H-bond dashes
cmd.distance("contacts", "ligand", "binding_site", cutoff=3.5)
cmd.show("dashes", "contacts")
cmd.color("forest", "contacts")
cmd.hide("labels", "contacts")
cmd.set("dash_radius", 0.06)

# Residue labels with Arial — keep the 6 closest to ligand (clean journal style)
lig_model = cmd.get_model("ligand")
lig_coords = [a.coord for a in lig_model.atom]
lcx = sum(c[0] for c in lig_coords) / len(lig_coords)
lcy = sum(c[1] for c in lig_coords) / len(lig_coords)
lcz = sum(c[2] for c in lig_coords) / len(lig_coords)

bs_atoms = cmd.get_model("binding_site and name CA").atom
# Sort residues by distance from CA to ligand center; keep 6 nearest
bs_atoms_sorted = sorted(bs_atoms, key=lambda a: (a.coord[0]-lcx)**2 + (a.coord[1]-lcy)**2 + (a.coord[2]-lcz)**2)
bs_atoms_keep = bs_atoms_sorted[:6]

# Center of mass of kept residues for outward offset direction
cmx = sum(a.coord[0] for a in bs_atoms_keep) / len(bs_atoms_keep)
cmy = sum(a.coord[1] for a in bs_atoms_keep) / len(bs_atoms_keep)
cmz = sum(a.coord[2] for a in bs_atoms_keep) / len(bs_atoms_keep)

for a in bs_atoms_keep:
    resn, resi = a.resn, a.resi
    x, y, z = a.coord
    label_name = f"label_{resi}"
    ca_sel = f"binding_site and resi {resi} and name CA"
    
    dx, dy, dz = x - cmx, y - cmy, z - cmz
    norm = max(0.001, (dx**2 + dy**2 + dz**2)**0.5)
    label_offset = 6.0
    label_pos = [x + dx/norm*label_offset, y + dy/norm*label_offset + 1.5, z + dz/norm*label_offset]
    
    cmd.pseudoatom(label_name, pos=label_pos)
    cmd.label(label_name, f'"{resn} {resi}"')
    cmd.hide("nonbonded", label_name)
    
    line_name = f"line_{resi}"
    cmd.distance(line_name, label_name, ca_sel, cutoff=99.0)
    cmd.show("dashes", line_name)
    cmd.color("gray60", line_name)
    cmd.set("dash_radius", 0.04, line_name)
    cmd.hide("labels", line_name)

try:
    cmd.set("label_font_name", "Arial")
except:
    cmd.set("label_font_id", 7)
cmd.set("label_size", 18)
cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")

# Optimal view
cmd.center("ligand")
cmd.zoom("ligand", 6)
cmd.turn("y", -35)
cmd.turn("x", 25)

# 3000 DPI rendering
cmd.ray(3600, 2800)
panel_b = os.path.join(OUT_DIR, "pravastatin_panelB_3000dpi.png")
cmd.png(panel_b, dpi=3000, width=3600, height=2800)
print(f"Panel B done: {panel_b}")

cmd.quit()
print("=== Pravastatin 3000 DPI render complete ===")
