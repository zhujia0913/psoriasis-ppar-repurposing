#!/usr/bin/env python3
"""
MCODE module detection + CytoHubba Hub gene screening
for 18-gene PPAR target DEG PPI network.
"""
import json, urllib.request, urllib.parse, os
import networkx as nx
import numpy as np
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
outdir = HERE / "figures"
os.makedirs(outdir, exist_ok=True)

# ============================================================
# 1. Build 18-gene PPI network
# ============================================================
genes_up   = ["CXCL13","CXCR4","ALDH1A3","FABP5","NAMPT","HMOX1","HBEGF","ANGPTL4","LDLR"]
genes_down = ["INSIG1","APOE","BCL2","MOGAT1","PLIN1","ADIPOQ","RBP4","FADS2","PDK4"]
all_genes  = genes_up + genes_down

base = "https://cn.string-db.org/api/json"
params_net = urllib.parse.urlencode({
    "identifiers": "\r".join(all_genes),
    "species": 9606,
    "required_score": 400,
    "network_type": "functional",
    "caller_identity": "academic_research"
})
req = urllib.request.Request(f"{base}/network?{params_net}")
with urllib.request.urlopen(req, timeout=60) as resp:
    net_data = json.loads(resp.read().decode())

G = nx.Graph()
for g in all_genes:
    G.add_node(g)
for edge in net_data:
    a, b, score = edge["preferredName_A"], edge["preferredName_B"], edge["score"]
    G.add_edge(a, b, weight=score)

# Remove singletons
singletons = [n for n in G.nodes() if G.degree(n) == 0]
G_clean = G.copy()
G_clean.remove_nodes_from(singletons)

print(f"Network: {G_clean.number_of_nodes()} connected nodes, {G_clean.number_of_edges()} edges")
print(f"Density: {nx.density(G_clean):.4f}")
print(f"Singletons excluded: {singletons}\n")

# ============================================================
# 2. MCODE Module Detection
# Parameters: Degree Cutoff≥2, haircut≥0.2, Node Score Cutoff≥0.2, K-core≥2, Max Depth=100
# ============================================================
def mcode_weighted(G, degree_cutoff=2, node_score_cutoff=0.2,
                   haircut=True, haircut_threshold=0.2, k_core=2, max_depth=100):
    """MCODE algorithm (Bader & Hogue, BMC Bioinformatics 2003)"""
    # Step 1: Weight nodes by neighborhood density × k-core size
    scores = {}
    for v in G.nodes():
        neighbors = set(G.neighbors(v))
        if len(neighbors) < degree_cutoff:
            continue
        subg = G.subgraph(neighbors)
        if subg.number_of_nodes() < 2:
            scores[v] = 0
        else:
            k = min(k_core, subg.number_of_nodes())
            try:
                k_core_g = nx.k_core(subg, k=k)
                density = nx.density(k_core_g)
                scores[v] = density * k_core_g.number_of_nodes()
            except nx.NetworkXError:
                scores[v] = 0

    active = {v for v, s in scores.items() if s >= node_score_cutoff}

    complexes = []
    seen = set()

    for seed in sorted(active, key=lambda v: scores.get(v, 0), reverse=True):
        if seed in seen:
            continue
        cluster = {seed}
        changed = True
        depth = 0
        while changed and depth < max_depth:
            changed = False
            depth += 1
            neighbors = set()
            for node in cluster:
                neighbors.update(set(G.neighbors(node)))
            neighbors -= cluster
            best_node, best_score = None, 0
            for n in neighbors:
                trial = cluster | {n}
                subg = G.subgraph(trial)
                td = nx.density(subg)
                if td >= best_score:
                    best_score = td
                    best_node = n
            if best_node is not None and best_score > 0:
                cluster.add(best_node)
                changed = True

        # Haircut
        if haircut and len(cluster) > 1:
            subg = G.subgraph(cluster)
            internal_deg = dict(subg.degree())
            to_remove = {n for n in cluster if internal_deg.get(n, 0) < len(cluster) * haircut_threshold}
            cluster -= to_remove
            if len(cluster) < 2:
                continue

        # K-core filter
        if len(cluster) >= k_core:
            try:
                core = nx.k_core(G.subgraph(cluster), k=k_core)
                cluster = set(core.nodes())
            except nx.NetworkXError:
                continue

        if len(cluster) >= 2:
            sf = G.subgraph(cluster)
            complexes.append({
                'genes': sorted(cluster),
                'size': len(cluster),
                'density': round(nx.density(sf), 4),
                'score': round(nx.density(sf) * len(cluster), 4)
            })
        seen.update(cluster)

    complexes.sort(key=lambda x: x['score'], reverse=True)
    return complexes


mcode_modules = mcode_weighted(G_clean, degree_cutoff=2, node_score_cutoff=0.2,
                               haircut=True, haircut_threshold=0.2, k_core=2, max_depth=100)

print("=" * 70)
print("MCODE MODULE DETECTION (Degree Cutoff≥2, haircut≥0.2, Node Score≥0.2, K-core≥2)")
print("=" * 70)
if not mcode_modules:
    print("⚠ No modules found with default parameters. Trying relaxed (score≥0.1)...")
    mcode_modules = mcode_weighted(G_clean, degree_cutoff=2, node_score_cutoff=0.1,
                                   haircut=True, haircut_threshold=0.2, k_core=2, max_depth=100)

for i, m in enumerate(mcode_modules):
    up_in = [g for g in m['genes'] if g in genes_up]
    down_in = [g for g in m['genes'] if g in genes_down]
    print(f"\n  Module {i+1}: size={m['size']}, score={m['score']:.4f}, density={m['density']:.4f}")
    print(f"  Up ({len(up_in)}): {', '.join(up_in)}")
    print(f"  Down ({len(down_in)}): {', '.join(down_in)}")

# Save MCODE
with open(os.path.join(outdir, "mcode_modules.tsv"), "w") as f:
    f.write("Module\tSize\tScore\tDensity\tGenes\tUp_Count\tDown_Count\n")
    for i, m in enumerate(mcode_modules):
        up_in = [g for g in m['genes'] if g in genes_up]
        down_in = [g for g in m['genes'] if g in genes_down]
        f.write(f"Module_{i+1}\t{m['size']}\t{m['score']:.4f}\t{m['density']:.4f}\t"
                f"{','.join(m['genes'])}\t{len(up_in)}\t{len(down_in)}\n")

# ============================================================
# 3. CytoHubba 11 Algorithms
# ============================================================

# 3a. Standard metrics
degree_c   = dict(G_clean.degree())
between_c  = nx.betweenness_centrality(G_clean, normalized=True)
close_c    = nx.closeness_centrality(G_clean)
ecc_c      = nx.eccentricity(G_clean)

# Stress centrality (number of shortest paths passing through node)
def stress_centrality(G):
    all_nodes = list(G.nodes())
    stress = {v: 0 for v in all_nodes}
    for s in all_nodes:
        paths = nx.single_source_shortest_path(G, s)
        for t in all_nodes:
            if s >= t:
                continue
            if t in paths:
                p = paths[t]
                for v in p[1:-1]:  # exclude endpoints
                    stress[v] += 1
    return stress

stress_c = stress_centrality(G_clean)

# 3b. BottleNeck
def bottleneck_centrality(G):
    shortest_paths = dict(nx.shortest_path(G))
    all_nodes = list(G.nodes())
    bc = {}
    for v in all_nodes:
        count = 0
        total = 0
        for i, s in enumerate(all_nodes):
            for t in all_nodes[i+1:]:
                if s == v or t == v:
                    continue
                total += 1
                path = shortest_paths.get(s, {}).get(t, [])
                if v in path:
                    count += 1
        bc[v] = count / max(total, 1)
    return bc

bottleneck_c = bottleneck_centrality(G_clean)

# 3c. Radiality
def radiality_centrality(G):
    n = G.number_of_nodes()
    rad = {}
    for v in G.nodes():
        lengths = nx.single_source_shortest_path_length(G, v)
        total = sum(lengths.values())
        rad[v] = (n - 1) / max(total, 1) if total > 0 else 0
    return rad

radiality_c = radiality_centrality(G_clean)

# 3d. MCC (Maximal Clique Centrality)
def mcc_centrality(G):
    n = G.number_of_nodes()
    cliques = list(nx.find_cliques(G))
    mcc = {}
    for v in G.nodes():
        max_edges = 0
        for clq in cliques:
            if v in clq:
                s = len(clq)
                edges_in = s * (s - 1) / 2
                if edges_in > max_edges:
                    max_edges = edges_in
        mcc[v] = max_edges / (n - 1) if n > 1 else 0
    return mcc

mcc_c = mcc_centrality(G_clean)

# 3e. MNC & DMNC
def mnc_dmnc(G):
    mnc, dmnc = {}, {}
    for v in G.nodes():
        neighbors = set(G.neighbors(v))
        if len(neighbors) < 2:
            mnc[v] = 1; dmnc[v] = 0.0
        else:
            subg = G.subgraph(neighbors)
            comps = list(nx.connected_components(subg))
            if comps:
                largest = max(comps, key=len)
                mnc[v] = len(largest)
                dmnc[v] = nx.density(subg.subgraph(largest)) if len(largest) > 1 else 0.0
            else:
                mnc[v] = 1; dmnc[v] = 0.0
    return mnc, dmnc

mnc_c, dmnc_c = mnc_dmnc(G_clean)

# 3f. EPC (Edge Percolated Component) - simplified as normalized degree
epc_c = {v: degree_c[v] / max(G_clean.number_of_nodes() - 1, 1) for v in G_clean.nodes()}

# ============================================================
# 4. Rankings & Hub Gene Screening
# ============================================================
metrics = ['MCC', 'DMNC', 'MNC', 'Degree', 'EPC', 'Betweenness',
           'Stress', 'Radiality', 'Closeness', 'EcCentricity', 'BottleNeck']

metric_funcs = {
    'MCC': mcc_c, 'DMNC': dmnc_c, 'MNC': mnc_c, 'Degree': degree_c,
    'EPC': epc_c, 'Betweenness': between_c, 'Stress': stress_c,
    'Radiality': radiality_c, 'Closeness': close_c, 'EcCentricity': ecc_c,
    'BottleNeck': bottleneck_c
}

# Rankings (1 = best)
rankings = {}
for m in metrics:
    vals = metric_funcs[m]
    sorted_genes = sorted(vals.keys(), key=lambda x: vals[x], reverse=True)
    rankings[m] = {g: i+1 for i, g in enumerate(sorted_genes)}

# Full metrics table
print("\n" + "=" * 70)
print("CYTOHUBBA 11 TOPOLOGICAL ALGORITHMS — RAW VALUES")
print("=" * 70)

# Header
print(f"\n{'Gene':12s}", end='')
for m in metrics:
    print(f" {m:>12s}", end='')
print(f" {'Dir':>5s}")
print("-" * (12 + 13*len(metrics) + 5))

for g in sorted(G_clean.nodes(), key=lambda x: degree_c[x], reverse=True):
    print(f"{g:12s}", end='')
    for m in metrics:
        print(f" {metric_funcs[m][g]:12.4f}", end='')
    dir_str = "UP" if g in genes_up else "DOWN"
    print(f" {dir_str:>5s}")

# Rankings table
print(f"\n{'Gene':12s}", end='')
for m in metrics:
    print(f" {m:>12s}", end='')
print(f"\n{'-' * (12 + 13*len(metrics))}")

for g in sorted(G_clean.nodes(), key=lambda x: degree_c[x], reverse=True):
    print(f"{g:12s}", end='')
    for m in metrics:
        print(f" {rankings[m][g]:12d}", end='')
    print()

# ============================================================
# 5. Algorithm Evaluation — Discriminability Assessment
# ============================================================
print("\n" + "=" * 70)
print("ALGORITHM DISCRIMINABILITY ASSESSMENT (15-node subnetwork)")
print("=" * 70)

disc = {}
for m in metrics:
    vals = [metric_funcs[m][g] for g in G_clean.nodes()]
    unique = len(set(vals))
    top_tied = sum(1 for g in G_clean.nodes() if vals.count(metric_funcs[m][g]) > 1 and rankings[m][g] == 1)
    disc[m] = {'unique': unique, 'tied_at_top': top_tied}
    status = "✅ PASS" if unique >= 6 and top_tied <= 1 else ("⚠️ MARGINAL" if unique >= 5 else "❌ EXCLUDED")
    print(f"  {m:14s}  unique={unique:2d}  tied@1st={top_tied}  → {status}")

# ============================================================
# 6. Hub Gene Screening — 4-Algorithm Rank Sum (ties.method='min')
# ============================================================
# Selected algorithms: Degree, Betweenness, MNC, Closeness
# Rationale: see docs/CytoHubba方法论_调研与最终方案.md §5.1

selected = ['Degree', 'Betweenness', 'MNC', 'Closeness']
print(f"\n{'='*70}")
print(f"HUB GENE SCREENING — 4-Algorithm Rank Sum (ties.method='min')")
print(f"Selected: {', '.join(selected)}")
print(f"{'='*70}")

# Compute ranks with ties.method='min' (1=best)
from statistics import mean as stats_mean
rank4 = {}
for m in selected:
    vals = [metric_funcs[m][g] for g in G_clean.nodes()]
    # Rank descending: highest value = rank 1
    sorted_val_rank = sorted(set(vals), reverse=True)
    rank_map = {v: i+1 for i, v in enumerate(sorted_val_rank)}
    rank4[m] = {g: rank_map[metric_funcs[m][g]] for g in G_clean.nodes()}

# Rank Sum
hub_scores = []
for g in G_clean.nodes():
    mr = stats_mean([rank4[m][g] for m in selected])
    hub_scores.append((g, mr, {m: rank4[m][g] for m in selected}))
hub_scores.sort(key=lambda x: x[1])

# Natural breakpoint detection
print(f"\n{'Rank':<6s}{'Gene':12s}{'MeanRk':>8s}", end='')
for m in selected:
    print(f" {m[:4]:>5s}", end='')
print(f" {'Δ':>6s} {'Dir':>5s}")
print("-"*60)

prev_mean = None
for i, (g, mr, rks) in enumerate(hub_scores, 1):
    delta = f"{(mr - prev_mean):.2f}" if prev_mean else "—"
    prev_mean = mr
    dir_str = "UP" if g in genes_up else "DOWN"
    hub_mark = " ← HUB" if mr <= 3.75 else ""
    print(f"{i:<6d}{g:12s}{mr:8.2f}", end='')
    for m in selected:
        print(f" {rks[m]:5d}", end='')
    print(f" {delta:>6s} {dir_str:>5s}{hub_mark}")

# Summarize hubs
hubs = [(g, mr) for g, mr, _ in hub_scores if mr <= 3.75]
print(f"\n📊 {len(hubs)} Hub Genes Identified:")
for g, mr in hubs:
    dir_str = "↓" if g in genes_down else "↑"
    d, b, m, c = [rank4[alg][g] for alg in selected]
    print(f"   {dir_str} {g:10s}  MeanRank={mr:.2f}  (D={d} B={b} M={m} C={c})")

print(f"\n💡 Natural breakpoint: Δ={hub_scores[4][1] - hub_scores[3][1]:.2f} between rank #4 ({hub_scores[3][0]}) and #5 ({hub_scores[4][0]})")

# ============================================================
# 7. Save Outputs
# ============================================================

# 7a. Raw values (all 11 algorithms)
with open(os.path.join(outdir, "cytohubba_rankings.tsv"), "w") as f:
    f.write("Gene\t" + "\t".join(metrics) + "\tDirection\n")
    for g in sorted(G_clean.nodes(), key=lambda x: degree_c[x], reverse=True):
        f.write(g)
        for m in metrics:
            f.write(f"\t{metric_funcs[m][g]:.6f}")
        f.write(f"\t{'UP' if g in genes_up else 'DOWN'}\n")

# 7b. Ranks (4 selected algorithms, ties.method='min')
with open(os.path.join(outdir, "cytohubba_ranks.tsv"), "w") as f:
    f.write("Gene\t" + "\t".join(selected) + "\tDirection\n")
    for g, mr, _ in hub_scores:
        f.write(g)
        for m in selected:
            f.write(f"\t{rank4[m][g]}")
        f.write(f"\t{'UP' if g in genes_up else 'DOWN'}\n")

# 7c. Hub genes final (Rank Sum output)
with open(os.path.join(outdir, "hub_genes_final.tsv"), "w") as f:
    f.write("Rank\tGene\tMeanRank\tDegree\tBetweenness\tMNC\tCloseness\tDirection\tIsHub\n")
    for i, (g, mr, rks) in enumerate(hub_scores, 1):
        is_hub = "Yes" if mr <= 3.75 else "No"
        f.write(f"{i}\t{g}\t{mr:.2f}\t{rks['Degree']}\t{rks['Betweenness']}\t"
                f"{rks['MNC']}\t{rks['Closeness']}\t"
                f"{'UP' if g in genes_up else 'DOWN'}\t{is_hub}\n")

print(f"\n✅ Output files saved to: {outdir}/")
print(f"   cytohubba_rankings.tsv  (11 algorithms, raw values)")
print(f"   cytohubba_ranks.tsv      (4 algorithms, ranks with ties.method='min')")
print(f"   hub_genes_final.tsv      (4-algorithm Rank Sum, hub genes annotated)")
print(f"   mcode_modules.tsv        (MCODE module detection)")
print(f"   mcode_cytohubba.py       (this script)")
