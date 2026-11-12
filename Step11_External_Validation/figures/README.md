# figures/

Runtime output directory. `scripts/02_external_validation.R` writes
`Fig16_external_validation.pdf` and `Fig16_external_validation.png` (5908 ×
5413 px, 600 dpi) here. These regenerated files are git-ignored; a
journal-ready 600-dpi LZW TIFF reference copy is committed at
`../results/Fig16_external_validation.tif`.

To reproduce the archived TIFF from the regenerated PNG:

```bash
python3 -c "from PIL import Image; im = Image.open('figures/Fig16_external_validation.png').convert('RGB'); im.save('results/Fig16_external_validation.tif', compression='tiff_lzw', dpi=(600,600))"
```

The script calls `set.seed(20261111)` before drawing so the jittered points
are deterministic and the TIFF can be regenerated pixel-for-pixel.
