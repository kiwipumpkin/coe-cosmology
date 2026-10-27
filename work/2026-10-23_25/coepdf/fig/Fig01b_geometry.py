# AUTO-WRAP generator for Fig01b_geometry.png
"""Fig01b_geometry.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: coepdf/scripts/fig01b_geometry.py (saves fig/Fig01b_geometry.png)
Embedded source: /workspace/coepdf/scripts/fig01b_geometry.py  (verbatim below, unmodified)
Run:  python3 Fig01b_geometry.py   -> writes Fig01b_geometry_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'Fig01b_geometry' + '_regen.png')
SRC_CWD = '/workspace/coepdf'
MATCH = ['Fig01b_geometry.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_Fig01b_geometry_')
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
sys.argv = ['/workspace/coepdf/scripts/fig01b_geometry.py'] + []

SOURCE = "# Fig 1 bottom: observer O (Earth) inside thick shell, distance D from CoE core, vision sphere R.\nimport matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,numpy as np\nfrom matplotlib.patches import Circle,FancyArrowPatch\nplt.rcParams.update({'mathtext.fontset':'cm','font.family':'serif','font.size':15})\nfig=plt.figure(figsize=(5.34,4.12),dpi=200,facecolor='white');ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,481);ax.set_ylim(371,0);ax.axis('off')\nC=(38,185);O=(214,185);D=O[0]-C[0]\nax.add_patch(Circle(C,262,fc='lightgrey',ec='0.35',lw=1.5))\nax.add_patch(Circle(C,96,fc='white',ec='0.35',lw=1.5))\nax.add_patch(Circle(C,D,fill=False,ec='tab:blue',ls='--',lw=2))\nax.add_patch(Circle(O,76,fc='tab:green',alpha=.3,ec='tab:green',lw=1.5))\nar=lambda a,b,c,lw=1.6:ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,color=c,lw=lw))\n# D: CoE -> observer, drawn as chord-free radius at an angle ending on dashed circle (|D| = CoE-to-O distance)\nar((C[0]+6,C[1]),(O[0]-7,O[1]),'tab:blue',2);ax.text(118,C[1]-8,'$D$',color='tab:blue',fontsize=18,ha='center')\nt2=np.radians(125);ar(O,(O[0]+76*np.cos(t2),O[1]-76*np.sin(t2)),'tab:green',2)\nax.text(O[0]-42,O[1]-20,'$R$',color='tab:green',fontsize=17)\nar(O,(472,185),'#0055cc',2.5);ax.text(390,218,r'$V_O\\!\\approx\\!7.7c$',color='#0055cc',fontsize=15)\nax.plot(*C,'o',color='tab:red',ms=10);ax.plot(*O,'o',color='#0055cc',mec='k',ms=9)\nax.text(8,C[1]+30,'CoE core',color='tab:red',ha='left',fontsize=15)\nax.text(68,118,'empty\\ninterior',color='0.4',ha='center',fontsize=14,linespacing=1.1)\nax.text(O[0]+6,O[1]+40,'Earth\\n(observer, $O$)',color='k',fontsize=13,ha='left',linespacing=1.05,bbox=dict(fc='white',ec='none',alpha=.75,pad=1))\nax.text(200,40,'shell',rotation=-42,color='0.3',fontsize=15,ha='center',va='center')\nax.text(270,40,'core $r_s$ = 163 Gpc\\n(4.9 D), off scale',color='tab:red',fontsize=14)\nax.add_patch(plt.Rectangle((0,0),481,371,fill=False,ec='k',lw=2))\nfig.savefig('fig/Fig01b_geometry.png',facecolor='white')\n"
exec(compile(SOURCE, '/workspace/coepdf/scripts/fig01b_geometry.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/scripts/fig01b_geometry.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
