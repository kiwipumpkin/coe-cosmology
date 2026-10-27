# AUTO-WRAP generator for recon_n3_2x2_pre_modes.png
"""recon_n3_2x2_pre_modes.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: pre-"modes" version of recon_n3_2x2.py saved as recon_n3_2x2.py.bak_modes
Embedded source: /workspace/cosmology/figures/fourier_resid/recon_n3_2x2.py.bak_modes  (verbatim below, unmodified)
Run:  python3 recon_n3_2x2_pre_modes.py   -> writes recon_n3_2x2_pre_modes_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/cosmology/figures/fourier_resid so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'recon_n3_2x2_pre_modes' + '_regen.png')
SRC_CWD = '/workspace/cosmology/figures/fourier_resid'
MATCH = ['recon_n3_2x2.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_recon_n3_2x2_pre_modes_')
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
sys.argv = ['/workspace/cosmology/figures/fourier_resid/recon_n3_2x2.py.bak_modes'] + []

SOURCE = '# Fig 25 2x2: top = data + CoE before + after n=1-3; bottom = residuals after +n=1, +n=1-2, n=1-3 overlaid.\nimport sys;sys.path.insert(0,\'..\')\nexec(open(\'recon_n3.py\').read().split("cols=[\'purple\'")[0])\nfrom matplotlib.lines import Line2D\nplt.rcParams.update({\'font.size\':14})\nRC=[\'tab:blue\',\'tab:orange\',\'magenta\'];lab=[\'after +n=1\',\'after +n=1–2\',\'after n=1–3\']\nfig=plt.figure(figsize=(16,9.6),facecolor=\'white\')\ngs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=0.3,wspace=0.16,left=.07,right=.99,top=.95,bottom=.16)\nbox=dict(fc=\'white\',ec=\'none\',alpha=.8)\ndef bands(x):x.axvspan(15,80,color=\'gold\',alpha=.18,zorder=0);x.axvspan(260,335,color=\'lightsteelblue\',alpha=.35,zorder=0)\nfor k,(nm,Y,sY,mb,mg) in enumerate(((\'A\',A0,sA,AFl,AFg),(\'B\',B0,sB,BLl,BLg))):\n  ms=[(mb.copy(),mg.copy())]\n  for n,p in steps2:ms.append((ms[-1][0]+T(n,p,th,k),ms[-1][1]+T(n,p,gt,k)))\n  st=lambda m:(np.sum((Y-m)**2/sY**2),np.corrcoef(Y,m)[0,1]**2)\n  a=fig.add_subplot(gs[0,k]);bands(a)\n  a.errorbar(l,Y,yerr=sY,fmt=\'o\',color=\'k\',ms=3,alpha=.5,lw=.8)\n  a.plot(g,ms[0][1],color=\'purple\',lw=2.5);a.plot(g,ms[3][1],color=\'tab:green\',lw=2.5)\n  (c0,r0),(c3,r3)=st(ms[0][0]),st(ms[3][0])\n  a.text(.02,.05,f\'before (CoE): χ²={c0:.2f}, r²={r0:.3f}\\nafter n=1–3: χ²={c3:.2f}, r²={r3:.3f}\',transform=a.transAxes,fontsize=12.5,bbox=box)\n  a.set_title(f\'{nm}: data, CoE model before and after n=1–3\',fontweight=\'bold\',fontsize=15)\n  a.set_ylabel(\'A/A$_{all}$\' if nm==\'A\' else \'B/B$_{all}$\');a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])\n  b=fig.add_subplot(gs[1,k]);bands(b);b.axhline(0,color=\'grey\',lw=1)\n  hh=[]\n  for j,i in enumerate((1,2,3)):\n    r=Y-ms[i][0];rms=np.sqrt(np.mean(r**2))\n    if i==3:b.errorbar(l,r,yerr=sY,fmt=\'none\',ecolor=\'0.6\',lw=.7,alpha=.6,zorder=1)\n    hh+=b.plot(l,r,\'-o\',color=RC[j],ms=3.5,lw=1.4,alpha=.85,label=f\'{lab[j]}: RMS={rms:.3f}\')\n  b.legend(handles=hh,loc=\'upper right\',fontsize=11.5,framealpha=.85)\n  b.set_title(f\'{nm} residuals (data − model), overlaid\',fontweight=\'bold\',fontsize=15)\n  b.set_ylabel(f\'{nm} − model\');b.set_xlim(0,360);b.set_xticks([0,90,180,270,360]);b.set_xlabel(\'ℓ [deg]\')\nh=[Line2D([],[],color=\'k\',marker=\'o\',ls=\'\',alpha=.5,label=\'reconstruction, 313 clusters (±σ Migkas Fig. 6; grey bars on n=1–3 residual)\'),Line2D([],[],color=\'purple\',lw=3,label=\'CoE model, before (A: AF_smooth, B: sheet BL)\'),Line2D([],[],color=\'tab:green\',lw=3,label=\'after shared n=1–3 removed\')]+[Line2D([],[],color=RC[j],marker=\'o\',lw=1.4,label=\'residual \'+lab[j]) for j in range(3)]+[plt.Rectangle((0,0),1,1,color=\'gold\',alpha=.3,label=\'ℓ 15–80°\'),plt.Rectangle((0,0),1,1,color=\'lightsteelblue\',alpha=.5,label=\'ℓ 260–335°\')]\nfig.legend(handles=h,loc=\'lower center\',ncol=3,frameon=False,fontsize=12.5)\nfig.savefig(\'recon_n3_2x2.png\',dpi=110,facecolor=\'white\')\n'
exec(compile(SOURCE, '/workspace/cosmology/figures/fourier_resid/recon_n3_2x2.py.bak_modes', 'exec'), {'__name__': '__main__', '__file__': '/workspace/cosmology/figures/fourier_resid/recon_n3_2x2.py.bak_modes'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
