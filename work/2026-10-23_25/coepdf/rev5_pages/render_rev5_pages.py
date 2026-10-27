"""Render script for rev5_pages/p-XX.png (24 pages): pdftoppm -r 60 -png of /workspace/coepdf/main_rev5.pdf.
(build.sh shows the same pdftoppm step at -r 45 for the contact sheet; these page PNGs are 510x660 = 60 dpi,
verified pixel-identical.) Writes to rev5_pages_regen/ beside this folder; never touches rev5_pages/.
Run: python3 render_rev5_pages.py"""
import os, subprocess, glob
H = os.path.dirname(os.path.abspath(__file__)); OUT = H + '_regen'; os.makedirs(OUT, exist_ok=True)
subprocess.run(['pdftoppm', '-r', '60', '-png', '/workspace/coepdf/main_rev5.pdf', os.path.join(OUT, 'p')], check=True)
print(len(glob.glob(OUT + '/p-*.png')), 'pages ->', OUT)
