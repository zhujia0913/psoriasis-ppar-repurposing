#!/usr/bin/env python3
"""可选工具：将 violin 位图转为 LZW 压缩 TIFF 并写入 300 DPI 元数据。

注意：plot_violin_literature_style.R 现已直接用 R 原生 tiff(LZW) 输出 300 DPI 的 TIFF，
本脚本仅作为手动重处理的备选工具，不再是主流程的一部分。

用法（在免疫浸润包根目录或 scripts/ 下运行）：
    python3 scripts/fix_dpi_tiff.py
"""
import os
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
src = os.path.join(FIG, "violin_cibersort_19cells_literature.png")
tiff_out = os.path.join(FIG, "violin_cibersort_19cells_literature.tiff")

if not os.path.exists(src):
    raise SystemExit(f"找不到源文件: {src}")

print("Opening PNG ...")
img = Image.open(src).convert("RGB")  # 白底合成，减小体积
print(f"  size: {img.size[0]}x{img.size[1]}, mode: {img.mode}")

# 300 DPI (期刊通用要求约 300 dpi；3000 dpi 会超出 10 MB 单文件上限)
print("Saving TIFF with LZW compression, 300 DPI ...")
img.save(tiff_out, format="TIFF", dpi=(300, 300), compression="tiff_lzw")
img.close()
print(f"✅ TIFF 已写入: {tiff_out}")
