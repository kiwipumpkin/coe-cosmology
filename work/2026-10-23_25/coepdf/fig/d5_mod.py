# AUTO-WRAP generator for d5_mod.png
"""d5_mod.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to coepdf/fig7/L16_d5_mod.png whose script is fig7/L16_d5_mod.py
Embedded source: /workspace/coepdf/fig7/L16_d5_mod.py  (verbatim below, unmodified)
Run:  python3 d5_mod.py   -> writes d5_mod_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'd5_mod' + '_regen.png')
SRC_CWD = '/workspace/coepdf'
MATCH = ['d5_mod.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_d5_mod_')
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
sys.argv = ['/workspace/coepdf/fig7/L16_d5_mod.py'] + []

SOURCE = "# PROVENANCE COPY of /workspace/coepdf/d5_mod_earth.py (original generating script, rev6 figure; run in its original directory).\n# Output figure copied here as L16_d5_mod.png. Not re-run for rev7: see Appendix B.\nimport numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt\nplt.rcParams.update({'font.size':11})\nexec(open('/workspace/coepdf/scripts/shp.py').read())\nr0,dv=42.0,14.3005\nm=lambda r:shp(r/42)\ns=(m(r0+1e-3)-m(r0-1e-3))/2e-3\nfig,ax=plt.subplots(figsize=(13,4.6));r=np.linspace(0,80,600)\nax.axvspan(r0-dv,r0+dv,color='gold',alpha=0.3,label=f'observable range $r_\\\\oplus\\\\pm${dv:.1f} Gpc ({r0-dv:.1f}–{r0+dv:.1f} Gpc)')\nre=r0+dv;mf=lambda r:np.where(r>re,m(re),m(r))\nax.axvspan(re,80,color='0.85',alpha=0.5,label='beyond observable edge: unconstrained')\nse=(m(re)-m(re-1e-3))/1e-3;x=np.linspace(0,80-re,200);ys=[]\nfor i,k in enumerate([2e-5,6e-5,1.5e-4,3e-4]):\n  y=m(re)+se*x-k*x**2;ys+=list(y);ax.plot(re+x,y,ls='--',lw=1.5,color=plt.cm.Blues(0.45+0.15*i),label='illustrative continuations (unconstrained)' if i==0 else None)\nri=r[r<=re];ax.plot(ri,m(ri),'b',lw=2,label=r'$m(r)=\\exp\\sum_k c_k(\\min(r/42,R_b)-1)^{k+1}$, $c$=[−0.041, −0.035, 0.1353, 0], $R_b$=2.011')\nrt=np.linspace(r0-1.6*dv,r0+1.6*dv,10);ax.plot(rt,m(r0)+s*(rt-r0),'r--',lw=1.6,label=f'tangent at Earth: slope {s:.2e} per Gpc')\nax.plot(r0,m(r0),'ko',ms=9,label=f'Earth, r = {r0:.0f} Gpc, m = {m(r0):.3f}')\nax.set_xlim(0,80);ax.set_ylim(0.97*min(min(m(ri)),min(ys)),1.03*max(max(m(ri)),max(ys)));ax.set_xlabel('distance from CoE axis point [Gpc]');ax.set_ylabel('m (before normalising)');ax.set_facecolor('w')\nfig.legend(loc='lower center',ncol=2,frameon=False,fontsize=9);fig.tight_layout(rect=(0,0.16,1,1));fig.savefig('fig/d5_mod.png',dpi=150,facecolor='w')\nprint(s,m(r0),m(r0-dv),m(r0+dv))\n"
exec(compile(SOURCE, '/workspace/coepdf/fig7/L16_d5_mod.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/fig7/L16_d5_mod.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
