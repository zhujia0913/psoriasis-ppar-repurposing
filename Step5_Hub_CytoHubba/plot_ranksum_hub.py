#!/usr/bin/env python3
"""
Hub基因Rank Sum可视化图
4算法(Degree, Betweenness, MNC, Closeness) Rank Sum + 自然断点分类
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os
from pathlib import Path

# === 读取数据 ===
HERE = Path(__file__).resolve().parent
data_path = HERE / "hub_genes_final.tsv"
df = pd.read_csv(data_path, sep='\t')

# 按MeanRank排序
df = df.sort_values('MeanRank').reset_index(drop=True)

# === 颜色映射 ===
colors = []
for _, row in df.iterrows():
    if row['IsHub'] == 'Yes':
        colors.append('#2196F3')  # Hub基因 - 蓝色
    else:
        colors.append('#BDBDBD')  # 非Hub - 灰色

# === Plot ===
fig, ax = plt.subplots(figsize=(10, 7), dpi=300)

# 水平柱状图（从下到上，rank越小越靠上）
y_pos = np.arange(len(df))
bars = ax.barh(y_pos, df['MeanRank'], color=colors, edgecolor='white', linewidth=0.8, height=0.7)

# 自然断点阈值线
threshold = 4.5  # BCL2(3.75)与CXCR4(5.75)之间
ax.axvline(x=threshold, color='#F44336', linestyle='--', linewidth=1.5, alpha=0.8, zorder=5)


# Y轴：基因名
ax.set_yticks(y_pos)
ax.set_yticklabels(df['Gene'].tolist(), fontsize=9)
# Hub基因的y轴标签加粗蓝色
hub_df = df[df['IsHub'] == 'Yes']
for i, (_, row) in enumerate(df.iterrows()):
    if row['IsHub'] == 'Yes':
        ax.get_yticklabels()[i].set_fontweight('bold')
        ax.get_yticklabels()[i].set_color('#1565C0')
    else:
        ax.get_yticklabels()[i].set_color('#9E9E9E')

# 柱子右侧标注4算法排名（仅Hub基因）
for i, (_, row) in enumerate(df.iterrows()):
    if row['IsHub'] == 'Yes':
        rank_detail = f"D={int(row['Degree'])}  B={int(row['Betweenness'])}  M={int(row['MNC'])}  C={int(row['Closeness'])}"
        ax.text(row['MeanRank'] + 0.1, i,
                f"  {rank_detail}",
                fontsize=7, va='center', color='#757575')

# X轴
ax.set_xlabel('Mean Rank (lower = more central)', fontsize=11)
ax.set_xlim(0, 15)
ax.invert_yaxis()  # rank 1 在顶部

# 标题
ax.set_title('Hub Gene Identification by CytoHubba Rank Sum\n(4 Algorithms: Degree + Betweenness + MNC + Closeness)', 
             fontsize=13, fontweight='bold', pad=15)

# 图例
hub_patch = mpatches.Patch(color='#2196F3', label='Hub gene (IsHub = Yes)')
nonhub_patch = mpatches.Patch(color='#BDBDBD', label='Non-hub')
threshold_line = plt.Line2D([0], [0], color='#F44336', linestyle='--', linewidth=1.5, label='Natural break threshold')
ax.legend(handles=[hub_patch, nonhub_patch, threshold_line], 
          loc='upper right', fontsize=9, framealpha=0.9)

# 美观
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

plt.tight_layout()

# === 保存 ===
output_path = HERE / "Figure10_RankSum_HubGenes.png"
plt.savefig(output_path, dpi=500, bbox_inches='tight', facecolor='white')
plt.close()

print(f"✅ 图已保存: {output_path}")
print(f"   Hub基因: {list(hub_df['Gene'])}")
print(f"   自然断点阈值: {threshold}")
