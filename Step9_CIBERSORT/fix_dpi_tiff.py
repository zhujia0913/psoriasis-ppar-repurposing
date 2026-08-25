#!/usr/bin/env python3
"""Fix PNG DPI metadata and create LZW-compressed TIFF for journal submission."""
from pathlib import Path
from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # disable decompression bomb check

HERE = Path(__file__).resolve().parent
figures = HERE / "figures"
src = str(figures / "violin_cibersort_22cells_literature.png")
tiff_out = str(figures / "violin_cibersort_22cells_literature.tiff")

print("Opening PNG ...")
img = Image.open(src)
print(f"  size: {img.size[0]}x{img.size[1]}, mode: {img.mode}")

# 1. Re-save PNG with 3000 DPI metadata
print("Saving PNG with 3000 DPI metadata ...")
img.save(src, dpi=(3000, 3000), optimize=False)

# 2. Save TIFF with LZW compression + 3000 DPI
print("Saving TIFF with LZW compression, 3000 DPI ...")
img.save(tiff_out, format="TIFF", dpi=(3000, 3000), compression="tiff_lzw")

img.close()
print("✅ All done")
