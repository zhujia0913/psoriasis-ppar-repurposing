#!/usr/bin/env python3
"""Fix PNG DPI metadata and create LZW-compressed TIFF for journal submission."""
import sys
from pathlib import Path
from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # disable decompression bomb check (our own 42000x36000 image)

HERE = Path(__file__).resolve().parent
figures = HERE / "figures"
src = str(figures / "Figure_GO_Enrichment_3x2.png")
png_out = src  # overwrite with correct DPI metadata
tiff_out = str(figures / "Figure_GO_Enrichment_3x2.tiff")

print("Opening PNG (42000x36000, ~4.5 GB in memory)...")
img = Image.open(src)
print(f"  size: {img.size[0]}x{img.size[1]}, mode: {img.mode}")

# 1. Re-save PNG with 3000 DPI metadata
print("Saving PNG with 3000 DPI metadata...")
img.save(png_out, dpi=(3000, 3000), optimize=False)
print("  PNG done")

# 2. Save as TIFF with LZW compression + 3000 DPI
print("Saving TIFF with LZW compression, 3000 DPI...")
img.save(tiff_out, format="TIFF", dpi=(3000, 3000), compression="tiff_lzw")
print("  TIFF done")

img.close()
print("✅ All done")
