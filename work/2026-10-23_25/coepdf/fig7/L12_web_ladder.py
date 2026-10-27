# AUTO-WRAP generator for L12_web_ladder.png
"""L12_web_ladder.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: restack of fig/f19mc3b.png (mc_br_plot MC web #3)
Embedded source: /workspace/coepdf/book/remade/f19mc3b_stack_remade.py  (verbatim below, unmodified)
Run:  python3 L12_web_ladder.py   -> writes L12_web_ladder_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/book/remade so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'L12_web_ladder' + '_regen.png')
SRC_CWD = '/workspace/coepdf/book/remade'
MATCH = ['f19mc3b_stack_remade.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_L12_web_ladder_')
import matplotlib; matplotlib.use('Agg')
import matplotlib.figure as _mf
_orig_open = builtins.open
def _divert(p):
    p = os.fspath(p)
    if os.path.abspath(p) == OUT: return p
    return os.path.join(SCRATCH, os.path.basename(p))
def _open(file, mode='r', *a, **k):
    if isinstance(file, (str, bytes, os.PathLike)) and any(c in mode for c in 'wax+'):
        file = _divert(file)
    return _orig_open(file, mode, *a, **k)
builtins.open = _open
_orig_save = _mf.Figure.savefig
def _savefig(self, fname, *a, **k):
    if isinstance(fname, (str, os.PathLike)):
        b = os.path.basename(os.fspath(fname))
        fname = OUT if b in MATCH else os.path.join(SCRATCH, b)
    return _orig_save(self, fname, *a, **k)
_mf.Figure.savefig = _savefig
try:
    from PIL import Image as _I
    _ois = _I.Image.save
    def _isave(self, fp, *a, **k):
        if isinstance(fp, (str, os.PathLike)):
            b = os.path.basename(os.fspath(fp)); fp = OUT if (b in MATCH or os.path.abspath(os.fspath(fp)) == OUT) else os.path.join(SCRATCH, b)
        return _ois(self, fp, *a, **k)
    _I.Image.save = _isave
except Exception: pass
os.chdir(SRC_CWD); sys.path.insert(0, SRC_CWD)
sys.argv = ['/workspace/coepdf/book/remade/f19mc3b_stack_remade.py'] + []

SOURCE = "# REMADE (not original): restack the 3x2 grid f19mc3b.png (from mc_br_plot.py, MC web #3) into one column.\nfrom PIL import Image\ng=Image.open('/workspace/coepdf/fig/f19mc3b.png').convert('RGB');W,H=g.size;w,h=W//3,H//2\ntiles=[g.crop((c*w,r*h,(c+1)*w,(r+1)*h)) for r in range(2) for c in range(3)]\nout=Image.new('RGB',(w,h*6),'white')\nfor i,t in enumerate(tiles):out.paste(t,(0,i*h))\nout.save('f19mc3b_stack_remade.png')\n"
exec(compile(SOURCE, '/workspace/coepdf/book/remade/f19mc3b_stack_remade.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/book/remade/f19mc3b_stack_remade.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
