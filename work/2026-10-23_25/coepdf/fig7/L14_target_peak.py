# AUTO-WRAP generator for L14_target_peak.png
"""L14_target_peak.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: book/remade/superseded/f13_remade.py (Unicode-label remake of Chris f13 code; also book/code/21_f13_remade.py)
Embedded source: /workspace/coepdf/book/remade/superseded/f13_remade.py  (verbatim below, unmodified)
Run:  python3 L14_target_peak.py   -> writes L14_target_peak_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/book/remade/superseded so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'L14_target_peak' + '_regen.png')
SRC_CWD = '/workspace/coepdf/book/remade/superseded'
MATCH = ['f13_remade.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_L14_target_peak_')
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
sys.argv = ['/workspace/coepdf/book/remade/superseded/f13_remade.py'] + []

SOURCE = "# REMADE (not original): non-axisymmetric source shell r=1.012R, Gaussian peak at l=220 (sigma 32, 8)\nimport numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt\nrng=np.random.default_rng(1);ell=np.arange(100,331)\nfig=plt.figure(figsize=(12,10));fig.suptitle('Non-axisymmetric source shell: density varies across latitude AND longitude\\nExterior shell r=1.012R • gravitational potential, not CMB temperature',fontsize=13)\nfor k,s in enumerate((32,8)):\n  ax=fig.add_subplot(2,2,2*k+1,projection='3d');u,v=np.mgrid[0:np.pi:120j,0:2*np.pi:240j]\n  f=rng.standard_normal(u.shape);clip=np.percentile(abs(f),97);col=plt.cm.RdBu_r((np.clip(f,-clip,clip)/clip+1)/2)\n  ax.plot_surface(np.sin(u)*np.cos(v),np.sin(u)*np.sin(v),np.cos(u),facecolors=col,rstride=1,cstride=1,linewidth=0,shade=False)\n  ax.set_box_aspect((1,1,1));ax.axis('off');ax.set_title(f'Distributed 3D clumps and voids\\nσℓ={s}; 97th-percentile color clip')\n  a=fig.add_subplot(2,2,2*k+2);t=np.exp(-.5*((ell-220)/s)**2)\n  a.plot(ell,t,'k',lw=2.2,label='Prescribed target');a.plot(ell,t,'--',color='#d24a3b',lw=1.8,label='Calculated spherical power')\n  a.axvline(220,color='.4',ls=':',lw=1);a.set(xlim=(110,330),ylim=(0,1.08),xlabel='Spherical multipole ℓ',ylabel='Normalized potential Dℓ',title=f'Peak ℓ=220, σℓ={s}');a.grid(alpha=.2);a.legend(fontsize=8)\nfig.tight_layout();fig.savefig('f13_remade.png',dpi=110,facecolor='white')\n"
exec(compile(SOURCE, '/workspace/coepdf/book/remade/superseded/f13_remade.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/book/remade/superseded/f13_remade.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
