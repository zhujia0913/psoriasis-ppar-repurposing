#!/usr/bin/env python3
"""
Step9 免疫浸润结果验证脚本 (独立复算, 不依赖 R) — BMC 投递版
用法: /usr/bin/python3 scripts/verify_step9_outputs.py
作用:
  1) 校验 CIBERSORT 细胞比例参数合法性 (形状/非负/每行和≈1)
  2) 用 scipy 独立重算 Hub×免疫 Spearman 相关 + BH-FDR, 与结果 tsv 逐项比对
  3) 检查是否存在 p==0 / padj==0 的数值下溢
  4) 报告关键对 (BCL2-Tregs)；若随包提供陈旧备份则量化差异
"""
import os, pandas as pd, numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # 免疫浸润/ 包根目录
PROP = os.path.join(BASE, "results", "cibersort_proportions_all.tsv")
TSV  = os.path.join(BASE, "results", "hub_immune_correlation_cibersort.tsv")
OLD  = os.path.join(BASE, "_PRE_RERUN_BACKUP_20261102", "results", "hub_immune_correlation_cibersort.tsv")  # 备份未随包提供
EXPR = os.path.join(BASE, "raw_data", "expression_pre_combat.csv")

print("===== A. CIBERSORT 比例参数合法性 =====")
prop = pd.read_csv(PROP, sep="\t", index_col=0)
noncell = [c for c in prop.columns if c.lower() in ("p-value","pvalue","correlation","pearson","rmse")]
cells = [c for c in prop.columns if c not in noncell]
frac = prop[cells]
print(f"  形状: {prop.shape}  (应为 209 样本 × 22 细胞)")
print(f"  含负数: {bool((frac<0).any().any())}   最小值: {frac.values.min():.4f}")
rs = frac.sum(axis=1)
print(f"  每行求和 min/mean/max: {rs.min():.4f}/{rs.mean():.4f}/{rs.max():.4f} (应≈1)")

print("\n===== B. 与独立 Python 复算逐项比对 (scipy + BH-FDR) =====")
expr = pd.read_csv(EXPR, index_col=0); expr.columns=[str(c) for c in expr.columns]
common=[s for s in prop.index if s in expr.columns]
hubs=["ADIPOQ","APOE","PLIN1","BCL2"]; cells2=list(cells)
rhos=[];ps=[]
for g in hubs:
    gv=expr.loc[g,common].to_numpy(float)
    for ct in cells2:
        r,p=spearmanr(gv,prop.loc[common,ct].to_numpy(float)); rhos.append(r);ps.append(p)
padj=multipletests(ps,method="fdr_bh")[1]
rec=pd.DataFrame({"Gene":np.repeat(hubs,len(cells2)),"CellType":np.tile(cells2,len(hubs)),
                  "Rho":rhos,"P.value":ps,"P.adj":padj})
old=pd.read_csv(TSV,sep="\t")
m=rec.merge(old,on=["Gene","CellType"],suffixes=("_r","_f"))
print(f"  rho 最大偏差: { (m['Rho_r']-m['Rho_f']).abs().max():.2e}")
print(f"  p   最大偏差: { (m['P.value_r']-m['P.value_f']).abs().max():.2e}")
print(f"  padj最大偏差: { (m['P.adj_r']-m['P.adj_f']).abs().max():.2e}")
print(f"  文件内 p==0 行数: {int((old['P.value']==0).sum())}   padj==0 行数: {int((old['P.adj']==0).sum())}")
bc=m[(m.Gene=="BCL2")&(m.CellType=="T cells regulatory (Tregs)")]
print(f"  BCL2-Tregs -> rho={bc['Rho_f'].iloc[0]:.4f} p={bc['P.value_f'].iloc[0]:.3e} padj={bc['P.adj_f'].iloc[0]:.3e}")

print("\n===== C. 相对陈旧备份的差异量化 (仅当备份随包提供时) =====")
if os.path.exists(OLD):
    oldv=pd.read_csv(OLD,sep="\t")
    mv=oldv.merge(old,on=["Gene","CellType"],suffixes=("_old","_new"))
    flip=int(((np.sign(mv['Rho_old'])!=0)&(np.sign(mv['Rho_new'])!=0)&(np.sign(mv['Rho_old'])!=np.sign(mv['Rho_new']))).sum())
    print(f"  配对总数: {len(mv)}")
    print(f"  符号翻转对数: {flip}")
    print(f"  |Δrho|>0.05 的对数: {int((mv['Rho_new']-mv['Rho_old']).abs().gt(0.05).sum())}")
    print(f"  旧版 P==0 行数: {int((mv['P.value_old']==0).sum())}  -> 新版: {int((mv['P.value_new']==0).sum())}")
else:
    print("  (陈旧备份未随包提供，跳过差异量化)")

print("\n===== 结论 =====")
ok = (m['Rho_r']-m['Rho_f']).abs().max()<1e-6 and (old['P.value']==0).sum()==0 and (old['P.adj']==0).sum()==0 and not (frac<0).any().any()
print("  PASS ✅" if ok else "  ⚠️ 需检查上述偏差")
