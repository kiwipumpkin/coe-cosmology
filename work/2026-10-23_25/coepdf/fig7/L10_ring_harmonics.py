# AUTO-WRAP generator for L10_ring_harmonics.png
"""L10_ring_harmonics.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: original f11 script never found; R10_harmonics.py is the computed remake
Embedded source: /workspace/coepdf/fig7/R10_harmonics.py  (verbatim below, unmodified)
Run:  python3 L10_ring_harmonics.py   -> writes L10_ring_harmonics_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/fig7 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'L10_ring_harmonics' + '_regen.png')
SRC_CWD = '/workspace/coepdf/fig7'
MATCH = ['R10_harmonics.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_L10_ring_harmonics_')
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
sys.argv = ['/workspace/coepdf/fig7/R10_harmonics.py'] + []

SOURCE = "# rev7 script for rev6 Fig 10 (axisymmetric density harmonics and their power). Background-independent geometry (no D or R enters).\n# Each panel: a ring of 400 bodies whose density sign follows P_l(cos θ) with 2% radial scatter; the bar chart is the power\n# actually measured from those bodies by Legendre projection (a_L = (2L+1)/2 ∫ s(θ) P_L sinθ dθ), not set by hand.\nimport numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt\nfrom scipy.special import eval_legendre\nrng=np.random.default_rng(0);th=np.linspace(0,2*np.pi,400,endpoint=False);mu=np.cos(th);t=np.linspace(0,np.pi,2001)\nfig=plt.figure(figsize=(9.7,13.7));fig.suptitle('Mass-density angular harmonics and their measured power\\nBlue: overdensity | Orange: underdensity',fontsize=13)\nfor L in range(8):\n    r,c=divmod(L,2);ax=fig.add_axes([0.03+c*0.5,0.765-r*0.24,0.27,0.19]);ax.set_title(f'ℓ = {L}, m = 0')\n    v=eval_legendre(L,mu);rad=1+0.02*rng.standard_normal(th.size)\n    ax.scatter(rad*np.sin(th),rad*np.cos(th),s=6,c=np.where(v>=0,'royalblue','darkorange'));ax.set_aspect('equal');ax.axis('off')\n    s=np.interp(t,np.sort(np.arccos(mu)),v[np.argsort(np.arccos(mu))])  # sampled density vs polar angle from the bodies\n    a=np.array([(2*k+1)/2*np.trapezoid(s*eval_legendre(k,np.cos(t))*np.sin(t),t) for k in range(8)]);p=a**2/(2*np.arange(8)+1);p/=p.max()\n    b=fig.add_axes([0.31+c*0.5,0.785-r*0.24,0.15,0.105]);b.bar(range(8),p,color='tab:blue')\n    b.set_xticks(range(8));b.set_yticks([0,1]);b.set_ylim(0,1.15);b.tick_params(labelsize=6);b.set_xlabel('Multipole L',fontsize=6)\nfig.savefig(__file__.replace('.py','.png'),dpi=110,facecolor='white')\n"
exec(compile(SOURCE, '/workspace/coepdf/fig7/R10_harmonics.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/fig7/R10_harmonics.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
