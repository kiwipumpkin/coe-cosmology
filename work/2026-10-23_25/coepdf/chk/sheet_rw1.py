"""sheet_rw1.png -- generating script. Original process: `coepdf/build.sh rw1` (tectonic main_rev5.tex ->
pdftoppm -r 45 -> 5-column PIL contact sheet saved as chk/sheet_rw1.png). The tex compiled at that moment is
saved as chk/main_rev5_rw1.tex. This script redoes those steps in a temp dir (build.sh itself would
rm -rf rev5_pages and overwrite main_rev5.pdf, so it is NOT run). Figures resolve via search path /workspace/coepdf.
Run: python3 sheet_rw1.py -> chk/sheet_rw1_regen.png"""
import os, subprocess, tempfile, shutil, glob
from PIL import Image
H = os.path.dirname(os.path.abspath(__file__)); T = tempfile.mkdtemp(prefix='sheet_rw1_')
shutil.copy(os.path.join(H, 'main_rev5_rw1.tex'), os.path.join(T, 'main_rev5.tex'))
subprocess.run(['tectonic', '-Z', 'search-path=/workspace/coepdf', 'main_rev5.tex'], cwd=T, check=True)
os.makedirs(T + '/rev5_pages'); subprocess.run(['pdftoppm', '-r', '45', '-png', T + '/main_rev5.pdf', T + '/rev5_pages/p'], check=True)
fs = sorted(glob.glob(T + '/rev5_pages/p-*.png')); ims = [Image.open(f) for f in fs]; w, h = ims[0].size
n = len(ims); c = 5; r = (n + c - 1) // c; S = Image.new('RGB', (w * c, h * r), 'white')
for i, im in enumerate(ims): S.paste(im, ((i % c) * w, (i // c) * h))
S.save(os.path.join(H, 'sheet_rw1_regen.png')); print('pages', n)
