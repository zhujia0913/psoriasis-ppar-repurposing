# figures/

Runtime output directory. `scripts/02_external_validation.R` writes
`Fig16_external_validation.pdf` and `Fig16_external_validation.png` (300 dpi)
here. These regenerated files are not committed; a journal-ready 300-dpi LZW
TIFF reference copy is available in `../results/Fig16_external_validation.tif`.

To convert the regenerated PNG to a TIFF:

```bash
python3 -c "from PIL import Image; im = Image.open('Fig16_external_validation.png'); im.save('Fig16_external_validation.tif', compression='tiff_lzw', dpi=(300,300))"
```
