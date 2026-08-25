#!/usr/bin/env python3
"""
PPI Network Visualization for MCODE Module 1 (10 genes).
Same style as the 18-gene PPI network:
- Hub genes marked with star shape and golden border
- Non-hub genes as circles colored by regulation direction
- Edge width proportional to STRING confidence score
- All labels in black below nodes
- Legend in upper right
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import networkx as nx
import numpy as np
import os
from pathlib import Path

# ============================================================
# Data
# ============================================================
HERE = Path(__file__).resolve().parent
outdir = HERE

# Module 1 genes
genes_up   = ["CXCL13", "CXCR4", "HMOX1", "ANGPTL4"]          # 4 up
genes_down = ["APOE", "BCL2", "PLIN1", "ADIPOQ", "RBP4", "PDK4"]  # 6 down

# Hub genes (4-algorithm Rank Sum top 4, natural breakpoint)
hub_genes = {"ADIPOQ", "APOE", "BCL2", "PLIN1"}

# logFC mapping
logfc = {
    "CXCL13": 2.66, "CXCR4": 1.60, "HMOX1": 1.19, "ANGPTL4": 1.02,
    "APOE": -1.04, "BCL2": -1.07, "PLIN1": -1.12, "ADIPOQ": -1.34,
    "RBP4": -1.36, "PDK4": -1.85
}

# Edges from STRING (score >= 0.4) — only within Module 1
edges = [
    ("PDK4", "PLIN1", 0.421), ("PDK4", "ADIPOQ", 0.440), ("PDK4", "ANGPTL4", 0.739),
    ("HMOX1", "ADIPOQ", 0.527), ("HMOX1", "APOE", 0.572), ("HMOX1", "BCL2", 0.856),
    ("APOE", "PLIN1", 0.403), ("APOE", "RBP4", 0.439), ("APOE", "ANGPTL4", 0.454),
    ("APOE", "BCL2", 0.491), ("APOE", "CXCL13", 0.557), ("APOE", "CXCR4", 0.558),
    ("APOE", "ADIPOQ", 0.795),
    ("CXCL13", "BCL2", 0.725), ("CXCL13", "CXCR4", 0.998),
    ("PLIN1", "ANGPTL4", 0.405), ("PLIN1", "ADIPOQ", 0.947),
    ("ANGPTL4", "ADIPOQ", 0.492),
    ("RBP4", "ADIPOQ", 0.799),
    ("BCL2", "ADIPOQ", 0.442), ("BCL2", "CXCR4", 0.759),
    ("CXCR4", "ADIPOQ", 0.527),
]

# ============================================================
# Build graph
# ============================================================
G = nx.Graph()
all_genes = genes_up + genes_down
for g in all_genes:
    G.add_node(g)

for a, b, s in edges:
    G.add_edge(a, b, weight=s)

singletons = [n for n in G.nodes() if G.degree(n) == 0]
print(f"Module 1: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
if singletons:
    print(f"Isolated nodes: {singletons}")
else:
    print("No isolated nodes (all connected)")

# ============================================================
# Layout
# ============================================================
pos = nx.spring_layout(G, k=2.5, iterations=300, seed=42, weight='weight')

# Manually adjust nodes that overlap with the legend (upper right area)
# Legend occupies roughly x>0.25, y>0.55
for n, (x, y) in list(pos.items()):
    if x > 0.2 and y > 0.5:
        pos[n] = (x, 0.45)  # push below legend area

# ============================================================
# Plot
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(10, 8), facecolor='white')
ax.set_facecolor('white')

# Colors
color_up = "#E74C3C"
color_down = "#3498DB"
hub_edge_color = "#F1C40F"
hub_marker = "*"
node_marker = "o"

# --- Draw edges ---
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
min_w, max_w = min(edge_weights), max(edge_weights)
edge_widths = [0.5 + 2.0 * (w - min_w) / (max_w - min_w) for w in edge_weights]

nx.draw_networkx_edges(G, pos, ax=ax,
                       edge_color='#BBBBBB',
                       width=edge_widths,
                       alpha=0.35)

# --- Draw non-hub nodes ---
non_hub_up = [n for n in G.nodes() if n not in hub_genes and n in genes_up]
non_hub_down = [n for n in G.nodes() if n not in hub_genes and n in genes_down]

for n in non_hub_up:
    x, y = pos[n]
    ax.scatter(x, y, s=1000, c=color_up, edgecolors='#333333',
              linewidths=0.8, zorder=3, marker=node_marker)
for n in non_hub_down:
    x, y = pos[n]
    ax.scatter(x, y, s=1000, c=color_down, edgecolors='#333333',
              linewidths=0.8, zorder=3, marker=node_marker)

# --- Draw hub nodes (star shape, golden border, larger) ---
hub_down = [n for n in G.nodes() if n in hub_genes and n in genes_down]
hub_up = [n for n in G.nodes() if n in hub_genes and n in genes_up]

for n in hub_down:
    x, y = pos[n]
    ax.scatter(x, y, s=1200, c=color_down, edgecolors=hub_edge_color,
              linewidths=2.5, zorder=4, marker=hub_marker)

for n in hub_up:
    x, y = pos[n]
    ax.scatter(x, y, s=1200, c=color_up, edgecolors=hub_edge_color,
              linewidths=2.5, zorder=4, marker=hub_marker)

# --- Draw labels (all black, below nodes) ---
for n in G.nodes():
    x, y = pos[n]
    fontsize = 8 if n in hub_genes else 7.5
    fontweight = 'bold' if n in hub_genes else 'normal'
    ax.text(x, y - 0.08, n, fontsize=fontsize, fontweight=fontweight,
            color='black', ha='center', va='top', zorder=5)

# --- Title ---
ax.set_title("PPI Network of MCODE Module 1\n(Module Score = 4.89, Density = 0.489, 10 nodes, 22 edges)",
             fontsize=13, fontweight='bold', pad=15)

# --- Legend ---
legend_elements = [
    mlines.Line2D([], [], color=color_up, marker='o', linestyle='None',
                  markersize=10, markeredgecolor='#333333', markeredgewidth=0.8,
                  label='Up-regulated (4 genes)'),
    mlines.Line2D([], [], color=color_down, marker='o', linestyle='None',
                  markersize=10, markeredgecolor='#333333', markeredgewidth=0.8,
                  label='Down-regulated (6 genes)'),
    mlines.Line2D([], [], color=color_down, marker='*', linestyle='None',
                  markersize=12, markeredgecolor=hub_edge_color, markeredgewidth=2.0,
                  label='Hub gene (4 genes)'),
    mlines.Line2D([], [], color='#999999', linewidth=1.5, alpha=0.5,
                  label='Interaction (STRING v12, score > 0.4)'),
]

ax.legend(handles=legend_elements, loc='upper right',
          fontsize=9, frameon=True, fancybox=True, shadow=False,
          edgecolor='#CCCCCC', framealpha=0.9)

# Remove axis
ax.axis('off')

plt.tight_layout()

# Save
out_png = os.path.join(outdir, "PPI_network_module1.png")
out_pdf = os.path.join(outdir, "PPI_network_module1.pdf")
fig.savefig(out_png, dpi=300, bbox_inches='tight', pad_inches=0.3, facecolor='white')
fig.savefig(out_pdf, bbox_inches='tight', pad_inches=0.3, facecolor='white')
plt.close()

print(f"\n✅ Saved:")
print(f"   {out_png}")
print(f"   {out_pdf}")

# Print hub gene summary
print(f"\n📊 Hub genes:")
for g in sorted(hub_genes, key=lambda x: -G.degree(x)):
    d = G.degree(g)
    print(f"   {g}: degree={d}, logFC={logfc[g]}, direction={'DOWN' if g in genes_down else 'UP'}")
