#!/usr/bin/env python3
"""TRRUST v2 TF Regulatory Network for 4 Hub Genes (ADIPOQ, APOE, PLIN1, BCL2).
Methodology: Chen et al. 2026, PLoS ONE (DOI: 10.1371/journal.pone.0338309)
Database: TRRUST v2 (Han et al. 2018, NAR), https://www.grnpedia.org/trrust/
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
from pathlib import Path

# ─── Configuration ───────────────────────────────────────────────
# All paths resolved relative to this script's own location.
HERE = Path(__file__).resolve().parent
STEP8_ROOT = HERE.parent
OUT_DIR = STEP8_ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DPI = 300
HUB_GENES = {"ADIPOQ", "APOE", "PLIN1", "BCL2"}

# ─── Load logFC from validated PPI nodes file ────────────────────
PPI_NODES_PATH = STEP8_ROOT / ".." / "Step4_PPI_MCODE" / "ppi_nodes.tsv"
if PPI_NODES_PATH.exists():
    df_ppi = pd.read_csv(PPI_NODES_PATH, sep="\t")
    logfc_map = dict(zip(df_ppi["Gene"], df_ppi["logFC"]))
    print(f"Loaded logFC from {PPI_NODES_PATH}: {len(logfc_map)} genes")
else:
    # Fallback (should not happen — PPI_nodes.tsv is the single source of truth)
    raise FileNotFoundError(f"PPI nodes file not found: {PPI_NODES_PATH}")

# ─── Load TRRUST data ────────────────────────────────────────────
df = pd.read_csv(str(STEP8_ROOT / "raw_data" / "trrust_rawdata_human.tsv"),
                 sep="\t", header=None,
                 names=["TF", "Target", "Mode", "PMID"])

df_hub = df[df["Target"].isin(HUB_GENES)].copy()
print(f"TRRUST raw interactions for 4 Hub genes: {len(df_hub)} rows")

# ─── Aggregate: collapse multiple modes per TF→Target pair ───────
agg_rows = []
for (tf, target), grp in df_hub.groupby(["TF", "Target"]):
    modes = grp["Mode"].unique()
    pmids = ";".join(grp["PMID"].unique())
    if len(modes) == 1:
        mode = modes[0]
    else:
        mode = "Mixed"
    n_pmids = len(set(pmid for p in grp["PMID"].str.split(";") for pmid in p if pmid))
    agg_rows.append({"TF": tf, "Target": target, "Mode": mode,
                      "n_PMIDs": n_pmids, "PMIDs": pmids})
df_agg = pd.DataFrame(agg_rows)
n_edges = len(df_agg)
n_tfs = df_agg["TF"].nunique()
n_raw = sum(df_hub.groupby(["TF", "Target"]).size())
print(f"Aggregated: {n_edges} edges, {n_tfs} TFs | Raw: {n_raw} interactions")
print(f"Mode distribution:\n{df_agg['Mode'].value_counts()}")

# ─── Build NetworkX graph ────────────────────────────────────────
G = nx.DiGraph()

# TF nodes
for tf in df_agg["TF"].unique():
    G.add_node(tf, type="TF", is_hub=False)

# Hub gene nodes — logFC from validated PPI source, NOT hardcoded
hub_info = {}
for gene in sorted(HUB_GENES):
    lfc = logfc_map.get(gene, None)
    if lfc is None:
        raise ValueError(f"logFC missing for {gene} in {PPI_NODES_PATH}")
    direction = "Up" if lfc > 0 else "Down"
    color = "#B2182B" if lfc > 0 else "#2166AC"
    hub_info[gene] = {"direction": direction, "logFC": round(lfc, 4), "color": color}
    print(f"  {gene}: logFC={lfc:.4f} ({direction})")

for gene, info in hub_info.items():
    G.add_node(gene, type="Hub", is_hub=True, **info)

# Edges
mode_colors = {"Activation": "#E41A1C", "Repression": "#377EB8",
               "Unknown": "#999999", "Mixed": "#984EA3"}
for _, row in df_agg.iterrows():
    G.add_edge(row["TF"], row["Target"], mode=row["Mode"],
               color=mode_colors.get(row["Mode"], "#999999"),
               n_pmids=row["n_PMIDs"])

total_nodes = G.number_of_nodes()
print(f"Graph: {total_nodes} nodes, {G.number_of_edges()} edges")

# ─── Layout: Circular (hub genes in center, TFs on outer ring) ───
hub_nodes_list = [n for n in G.nodes() if G.nodes[n].get("is_hub")]
tf_nodes_list = [n for n in G.nodes() if not G.nodes[n].get("is_hub")]

# Sort TFs by degree (descending) so high-degree TFs are evenly spaced
tf_nodes_sorted = sorted(tf_nodes_list, key=lambda n: G.degree(n), reverse=True)

pos = {}

# Hub genes: place in a tight cluster at center
n_hubs = len(hub_nodes_list)
hub_radius = 0.15
for i, gene in enumerate(sorted(hub_nodes_list)):
    angle = 2 * np.pi * i / n_hubs + np.pi / 2  # start from top
    pos[gene] = (hub_radius * np.cos(angle), hub_radius * np.sin(angle))

# TFs: place on outer circle, evenly spaced
n_tfs = len(tf_nodes_sorted)
tf_radius = 1.0
for i, tf in enumerate(tf_nodes_sorted):
    angle = 2 * np.pi * i / n_tfs + np.pi / 2  # start from top
    pos[tf] = (tf_radius * np.cos(angle), tf_radius * np.sin(angle))

# ─── Figure ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(20, 20))

for mode, color in mode_colors.items():
    edges = [(u, v) for u, v, d in G.edges(data=True) if d["mode"] == mode]
    if edges:
        widths = [0.8 + 0.3 * G[u][v].get("n_pmids", 1) for u, v in edges]
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=edges,
                               edge_color=color, alpha=0.5,
                               width=widths, arrows=True,
                               arrowsize=8, arrowstyle='-|>',
                               connectionstyle='arc3,rad=0.1')

# TF nodes
tf_sizes = [720 + 120 * G.degree(n) for n in tf_nodes_list]
nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=tf_nodes_list,
                       node_size=tf_sizes, node_color='#FDB863',
                       edgecolors='#E08214', linewidths=0.8,
                       node_shape='o', alpha=0.9)

# Hub nodes
hub_sizes = [4400 for _ in hub_nodes_list]
hub_colors = [G.nodes[n].get("color", "#2166AC") for n in hub_nodes_list]
nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=hub_nodes_list,
                       node_size=hub_sizes, node_color=hub_colors,
                       edgecolors='#333333', linewidths=2.5,
                       node_shape='o', alpha=1.0)

# TF labels — place outside the circle with slight offset
for n in tf_nodes_list:
    x, y = pos[n]
    # Radial offset: push label outward
    r = np.sqrt(x**2 + y**2)
    if r > 0:
        offset_x = x / r * 0.08
        offset_y = y / r * 0.08
    else:
        offset_x, offset_y = 0, 0
    ax.text(x + offset_x, y + offset_y, n,
            fontsize=16, fontweight='normal', color='#333333',
            ha='center', va='center', rotation=0)

# Hub labels (white on colored badge)
for n, (x, y) in pos.items():
    if n in hub_nodes_list:
        ax.text(x, y, n, fontsize=24, fontweight='bold', ha='center', va='center',
                color='white', bbox=dict(boxstyle='round,pad=0.3',
                facecolor=hub_info[n].get("color", "#2166AC"),
                edgecolor='white', alpha=0.95))

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='#E41A1C', lw=4, label='Activation'),
    Line2D([0], [0], color='#377EB8', lw=4, label='Repression'),
    Line2D([0], [0], color='#999999', lw=4, label='Unknown'),
    Line2D([0], [0], color='#984EA3', lw=4, label='Mixed'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FDB863',
           markeredgecolor='#E08214', markersize=20, label='TF (size ∝ degree)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2166AC',
           markeredgecolor='#333', markersize=28, label='Hub Gene (↓ Downregulated)'),
]
ax.legend(handles=legend_elements, loc='lower left', fontsize=18,
          framealpha=0.9, edgecolor='#cccccc')

# Title with computed counts (not hardcoded)
ax.set_title(f"TRRUST v2 Transcription Factor Regulatory Network\n"
             f"4 Hub Genes (ADIPOQ, APOE, PLIN1, BCL2) — {n_tfs} TFs, {n_edges} Regulatory Interactions",
             fontsize=30, fontweight='bold', pad=20)
ax.axis('off')
plt.tight_layout()

# Save
png_path = OUT_DIR / "TRRUST_TF_network_4hubs.png"
pdf_path = OUT_DIR / "TRRUST_TF_network_4hubs.pdf"
fig.savefig(png_path, dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(pdf_path, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"Saved: {png_path}")
print(f"Saved: {pdf_path}")
plt.close()

# ─── Save data files ─────────────────────────────────────────────
# Edge table
edges_out = []
for u, v, d in G.edges(data=True):
    edges_out.append({"TF": u, "Target": v, "Mode": d["mode"], "n_PMIDs": d["n_pmids"]})
pd.DataFrame(edges_out).sort_values(["Target", "TF"]).to_csv(
    OUT_DIR / "trrust_edges.tsv", sep="\t", index=False)

# Node table
nodes_out = []
for n in G.nodes():
    nd = G.nodes[n]
    nodes_out.append({
        "Node": n, "Type": nd["type"],
        "Degree": G.degree(n),
        "Direction": nd.get("direction", ""),
        "logFC": nd.get("logFC", ""),
    })
pd.DataFrame(nodes_out).sort_values(["Type", "Node"]).to_csv(
    OUT_DIR / "trrust_nodes.tsv", sep="\t", index=False)

# Full aggregated interaction table (with PMIDs)
df_agg_sorted = df_agg.sort_values(["Target", "TF"])
df_agg_sorted.to_csv(OUT_DIR / "trrust_interactions_aggregated.tsv", sep="\t", index=False)

# Per-hub summary
print("\n=== Per-Hub TF Summary ===")
for hub in sorted(HUB_GENES):
    tfs = df_agg[df_agg["Target"] == hub]["TF"].unique()
    n_act = len(df_agg[(df_agg["Target"] == hub) & (df_agg["Mode"] == "Activation")])
    n_rep = len(df_agg[(df_agg["Target"] == hub) & (df_agg["Mode"] == "Repression")])
    n_unk = len(df_agg[(df_agg["Target"] == hub) & (df_agg["Mode"] == "Unknown")])
    n_mix = len(df_agg[(df_agg["Target"] == hub) & (df_agg["Mode"] == "Mixed")])
    print(f"  {hub}: {len(tfs)} TFs (Act:{n_act} Rep:{n_rep} Unk:{n_unk} Mix:{n_mix})")

# Multi-hub TFs
print("\n=== Multi-Hub TFs ===")
tf_hubs = df_agg.groupby("TF")["Target"].apply(set).to_dict()
for tf, targets in sorted(tf_hubs.items()):
    if len(targets) > 1:
        print(f"  {tf} → {', '.join(sorted(targets))} ({len(targets)} hubs)")

print(f"\nAll outputs saved to: {OUT_DIR}")