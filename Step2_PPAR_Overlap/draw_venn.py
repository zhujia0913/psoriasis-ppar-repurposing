#!/usr/bin/env python3
"""
Venn diagram: PPARgene human (177) ∩ DEG (759) = 18
Output: figures/PPARgene_human_DEG_venn.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from pathlib import Path

# Data
ppargene = 177
deg = 759
overlap = 18

HERE = Path(__file__).resolve().parent
outdir = HERE / "figures"

# Figure
fig, ax = plt.subplots(figsize=(6, 5))
v = venn2(
    subsets=(ppargene - overlap, deg - overlap, overlap),
    set_labels=('PPAR gene', 'DEG'),
    set_colors=('#E64B35', '#4DBBD5'),
    alpha=0.7,
    ax=ax
)

# Style
for text in v.set_labels:
    if text: text.set_fontsize(12)
for text in v.subset_labels:
    if text: text.set_fontsize(13)

# Add a title
ax.set_title('PPARgene v2.0 Human ∩ Psoriasis DEG', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
outpath = os.path.join(outdir, 'PPARgene_human_DEG_venn.png')
plt.savefig(outpath, dpi=1000, bbox_inches='tight', facecolor='white')
plt.close()
print(f'✅ Saved: {outpath}')
