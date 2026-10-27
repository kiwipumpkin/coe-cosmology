"""L01a_geometry_top.png -- generating script.
Provenance: panel (a) of Fig01_geometry.png, which is written by
/workspace/cosmology/papers/spreadsheet_model/bubble/universe_thin.py (as ../universe_thin_drawing.png;
figs_letter_exact_python/makers/make_fig01_02_geometry_achart.py copies it to Fig01_geometry.png).
L01a_geometry_top.png is an exact pixel crop of that image: box x=6..383, y=41..250 (377x209),
found by exact sub-image search (mean abs diff 0). The crop was originally done by hand; it is scripted here.
Run: python3 L01a_geometry_top.py -> L01a_geometry_top_regen.png (original not overwritten).
universe_thin.py runs in its own dir; its file writes are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile, runpy
import matplotlib; matplotlib.use('Agg'); import matplotlib.figure as mf
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'L01a_geometry_top_regen.png')
BUBBLE = '/workspace/cosmology/papers/spreadsheet_model/bubble'; SCR = tempfile.mkdtemp(prefix='regen_L01a_')
_o = builtins.open
def _open(f, mode='r', *a, **k):
    if isinstance(f, (str, os.PathLike)) and any(c in mode for c in 'wax+'): f = os.path.join(SCR, os.path.basename(os.fspath(f)))
    return _o(f, mode, *a, **k)
builtins.open = _open
_s = mf.Figure.savefig
mf.Figure.savefig = lambda self, fn, *a, **k: _s(self, os.path.join(SCR, os.path.basename(os.fspath(fn))), *a, **k)
cwd = os.getcwd(); sys.path.insert(0, BUBBLE); os.chdir(BUBBLE)
try: runpy.run_path(os.path.join(BUBBLE, 'universe_thin.py'), run_name='__main__')
finally: os.chdir(cwd); builtins.open = _o
full = Image.open(os.path.join(SCR, 'universe_thin_drawing.png'))
full.crop((6, 41, 6 + 377, 41 + 209)).save(OUT); print('wrote', OUT)
