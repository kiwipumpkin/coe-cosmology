# AUTO-WRAP generator for recon_n3_final.png
"""recon_n3_final.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to cosmology/figures/fourier_resid/recon_n3_final.png, script recon_n3_final.py beside it
Embedded source: /workspace/cosmology/figures/fourier_resid/recon_n3_final.py  (verbatim below, unmodified)
Run:  python3 recon_n3_final.py   -> writes recon_n3_final_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/cosmology/figures/fourier_resid so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'recon_n3_final' + '_regen.png')
SRC_CWD = '/workspace/cosmology/figures/fourier_resid'
MATCH = ['recon_n3_final.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_recon_n3_final_')
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
sys.argv = ['/workspace/cosmology/figures/fourier_resid/recon_n3_final.py'] + []

SOURCE = '# Fig 25: per column (A | B): overlay of data + CoE (before) + after n=1-3; then residuals after +n=1, +n=1-2, n=1-3.\nimport sys;sys.path.insert(0,\'..\')\nexec(open(\'recon_n3.py\').read().split("cols=[\'purple\'")[0])\nfrom matplotlib.lines import Line2D\ncols=[\'purple\',\'tab:blue\',\'tab:red\',\'tab:green\']\nplt.rcParams.update({\'font.size\':14})\nfig=plt.figure(figsize=(16,12.7),facecolor=\'white\')\ngs=fig.add_gridspec(4,2,height_ratios=[1.5,1,1,1],hspace=0.45,wspace=0.16,left=.07,right=.99,top=.97,bottom=.105)\nbox=dict(fc=\'white\',ec=\'none\',alpha=.8)\ndef bands(x):x.axvspan(15,80,color=\'gold\',alpha=.18,zorder=0);x.axvspan(260,335,color=\'lightsteelblue\',alpha=.35,zorder=0)\nfor k,(nm,Y,sY,mb,mg) in enumerate(((\'A\',A0,sA,AFl,AFg),(\'B\',B0,sB,BLl,BLg))):\n  ms=[(mb.copy(),mg.copy())]\n  for n,p in steps2:ms.append((ms[-1][0]+T(n,p,th,k),ms[-1][1]+T(n,p,gt,k)))\n  st=lambda m:(np.sum((Y-m)**2/sY**2),np.corrcoef(Y,m)[0,1]**2)\n  a=fig.add_subplot(gs[0,k]);bands(a)\n  a.errorbar(l,Y,yerr=sY,fmt=\'o\',color=\'k\',ms=3,alpha=.5,lw=.8)\n  a.plot(g,ms[0][1],color=\'purple\',lw=2.5);a.plot(g,ms[3][1],color=\'tab:green\',lw=2.5)\n  (c0,r0),(c3,r3)=st(ms[0][0]),st(ms[3][0])\n  a.text(.02,.05,f\'before (CoE): χ²={c0:.2f}, r²={r0:.3f}\\nafter n=1–3: χ²={c3:.2f}, r²={r3:.3f}\',transform=a.transAxes,fontsize=12.5,bbox=box)\n  a.set_title(f\'{nm}: CoE model before and after n=1–3\',fontweight=\'bold\',fontsize=15)\n  a.set_ylabel(\'A/A$_{all}$\' if nm==\'A\' else \'B/B$_{all}$\');a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])\n  b0=None\n  for j,i in enumerate((1,2,3)):\n    b=fig.add_subplot(gs[j+1,k],sharey=b0);b0=b0 or b;bands(b)\n    r=Y-ms[i][0];rms=np.sqrt(np.mean(r**2))\n    b.errorbar(l,r,yerr=sY,fmt=\'o\',color=\'k\',ms=3,alpha=.5,lw=.8);b.axhline(0,color=\'grey\',lw=1)\n    if i<3:\n      n2,p2=steps2[i];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[i+1],lw=2.5)\n      t=f\'RMS={rms:.3f}; next: n={n2}, C={np.hypot(p2[2],p2[3]):.3f}\'\n    else:t=f\'RMS={rms:.3f}; final residual\'\n    b.text(.02,.06,t,transform=b.transAxes,fontsize=12,bbox=box)\n    b.set_title(f\'{nm} residual after \'+(\'+n=1\' if i==1 else \'+n=1–2\' if i==2 else \'n=1–3\'),fontweight=\'bold\',fontsize=14)\n    b.set_ylabel(f\'{nm} − model\');b.set_xlim(0,360);b.set_xticks([0,90,180,270,360])\n    if i==3:b.set_xlabel(\'ℓ [deg]\')\nh=[Line2D([],[],color=\'k\',marker=\'o\',ls=\'\',alpha=.5,label=\'reconstruction, 313 clusters (±σ Migkas Fig. 6)\'),Line2D([],[],color=\'purple\',lw=3,label=\'CoE model, before (A: AF_smooth, B: sheet BL)\'),Line2D([],[],color=\'tab:green\',lw=3,label=\'top: after n=1–3; residuals: next n=3\'),Line2D([],[],color=\'tab:red\',lw=3,label=\'next component n=2\'),plt.Rectangle((0,0),1,1,color=\'gold\',alpha=.3,label=\'ℓ 15–80°\'),plt.Rectangle((0,0),1,1,color=\'lightsteelblue\',alpha=.5,label=\'ℓ 260–335°\')]\nfig.legend(handles=h,loc=\'lower center\',ncol=3,frameon=False,fontsize=12.5)\nfig.savefig(\'recon_n3_final.png\',dpi=110,facecolor=\'white\')\n'
exec(compile(SOURCE, '/workspace/cosmology/figures/fourier_resid/recon_n3_final.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/cosmology/figures/fourier_resid/recon_n3_final.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
