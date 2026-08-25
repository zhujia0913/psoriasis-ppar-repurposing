#!/usr/bin/env python3
"""
Composite molecular docking figures — PLOS ONE compliant.

Figure 14: Fenofibrate-PPARα (Panel A + Panel B)
Figure 15: Pravastatin-PPARα (Panel A + Panel B)

PLOS ONE specs:
- TIFF, LZW compression, RGB, flattened
- Combination figure: 300 DPI
- Width 2.63-7.5 in -> use full 7.5 in (double column)
- Height <= 8.75 in
- Panel labels: Arial Bold, 12 pt at final physical size
- File < 10 MB
"""
from PIL import Image, ImageDraw, ImageFont
import os

from pathlib import Path
HERE = Path(__file__).resolve().parent
FIG_DIR = str(HERE / "figures")
OUT_DIR = str(HERE / "figures")       # PNG preview copy
PLOS_DIR = str(HERE / "plos_tiff")    # final PLOS ONE TIFFs
for d in (OUT_DIR, PLOS_DIR):
    os.makedirs(d, exist_ok=True)

DPI = 300
TARGET_W_IN = 7.5          # full double-column width
TARGET_W_PX = int(TARGET_W_IN * DPI)   # 2250 px

# Layout proportions (in inches at final size)
PAD_SIDE_IN = 0.10         # left/right margin
PAD_TOP_IN = 0.35          # room for panel labels
PAD_BOT_IN = 0.10
GAP_IN = 0.15              # gap between panels

LABEL_PT = 12              # PLOS: 8-12 pt
ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def composite_2panel_plos(name_a, name_b, out_name, label_a="A", label_b="B"):
    img_a = Image.open(os.path.join(FIG_DIR, name_a)).convert("RGB")
    img_b = Image.open(os.path.join(FIG_DIR, name_b)).convert("RGB")

    pad_side = int(PAD_SIDE_IN * DPI)
    pad_top = int(PAD_TOP_IN * DPI)
    pad_bot = int(PAD_BOT_IN * DPI)
    gap = int(GAP_IN * DPI)

    # each panel width
    panel_w = (TARGET_W_PX - 2 * pad_side - gap) // 2   # ~1070 px

    def fit(img):
        ratio = panel_w / img.width
        return img.resize((panel_w, int(round(img.height * ratio))), Image.LANCZOS)

    img_a = fit(img_a)
    img_b = fit(img_b)

    # equalize heights (crop-extend to max height, keep white bg)
    h = max(img_a.height, img_b.height)

    def pad_to(img, hh):
        if img.height == hh:
            return img
        out = Image.new("RGB", (img.width, hh), (255, 255, 255))
        out.paste(img, (0, (hh - img.height) // 2))
        return out

    img_a = pad_to(img_a, h)
    img_b = pad_to(img_b, h)

    total_h = pad_top + h + pad_bot
    canvas = Image.new("RGB", (TARGET_W_PX, total_h), (255, 255, 255))
    canvas.paste(img_a, (pad_side, pad_top))
    canvas.paste(img_b, (pad_side + panel_w + gap, pad_top))

    # panel labels: Arial Bold, 12 pt -> px = pt/72 * DPI
    draw = ImageDraw.Draw(canvas)
    font_px = int(round(LABEL_PT / 72.0 * DPI))
    try:
        font = ImageFont.truetype(ARIAL_BOLD, font_px)
    except Exception:
        font = ImageFont.load_default()

    x_a = pad_side
    x_b = pad_side + panel_w + gap
    y_lab = int(0.08 * DPI)   # ~0.08 in below top edge
    draw.text((x_a, y_lab), label_a, fill=(0, 0, 0), font=font)
    draw.text((x_b, y_lab), label_b, fill=(0, 0, 0), font=font)

    # save TIFF (LZW, 300 DPI) + PNG preview copy
    tiff_path = os.path.join(PLOS_DIR, out_name.replace(".png", ".tif"))
    canvas.save(tiff_path, "TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    png_path = os.path.join(OUT_DIR, out_name.replace("3000dpi", "plos"))
    canvas.save(png_path, "PNG", dpi=(DPI, DPI))

    w_in, h_in = TARGET_W_PX / DPI, total_h / DPI
    size_mb = os.path.getsize(tiff_path) / 1024 / 1024
    print(f"OK {out_name}")
    print(f"   TIFF: {tiff_path}")
    print(f"   {canvas.size[0]}x{canvas.size[1]} px  {w_in:.2f}x{h_in:.2f} in  "
          f"{DPI} DPI  LZW  {size_mb:.1f} MB")
    return tiff_path


print("=== Figure 14: Fenofibrate-PPARa (PLOS ONE) ===")
composite_2panel_plos(
    "fenofibrate_panelA_3000dpi.png",
    "fenofibrate_panelB_3000dpi.png",
    "Figure_14_Fenofibrate_PPARa_3000dpi.png",
    "A", "B",
)

print("\n=== Figure 15: Pravastatin-PPARa (PLOS ONE) ===")
composite_2panel_plos(
    "pravastatin_panelA_3000dpi.png",
    "pravastatin_panelB_3000dpi.png",
    "Figure_15_Pravastatin_PPARa_3000dpi.png",
    "A", "B",
)

print("\n=== PLOS composites complete ===")
