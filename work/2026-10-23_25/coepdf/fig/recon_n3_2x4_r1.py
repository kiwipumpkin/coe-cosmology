# AUTO-WRAP generator for recon_n3_2x4_r1.png
"""recon_n3_2x4_r1.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to cosmology/figures/fourier_resid/recon_n3_2x4_r1.png, script recon_n3_2x4_r1.py beside it
Embedded source: /workspace/cosmology/figures/fourier_resid/recon_n3_2x4_r1.py  (verbatim below, unmodified)
Run:  python3 recon_n3_2x4_r1.py   -> writes recon_n3_2x4_r1_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/cosmology/figures/fourier_resid so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'recon_n3_2x4_r1' + '_regen.png')
SRC_CWD = '/workspace/cosmology/figures/fourier_resid'
MATCH = ['recon_n3_2x4_r1.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_recon_n3_2x4_r1_')
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
sys.argv = ['/workspace/cosmology/figures/fourier_resid/recon_n3_2x4_r1.py'] + []

SOURCE = '# Term-by-term shared harmonic decomposition of A and B (chain: cum2.py -> recon_n3.py), with residuals beneath.\nimport sys;sys.path.insert(0,\'..\')\nexec(open(\'recon_n3.py\').read().split("cols=[\'purple\'")[0])\nfrom matplotlib.lines import Line2D\ncols=[\'purple\',\'tab:blue\',\'tab:red\',\'tab:green\']\nplt.rcParams.update({\'font.size\':19})\nfig=plt.figure(figsize=(26,17),facecolor=\'white\');gs=fig.add_gridspec(5,4,height_ratios=[1.25,0.55,0.16,1.25,0.55],hspace=0.08)\nAX={}\nstats={}\nfor blk,(nm,Y,sY,mb,mg,k) in enumerate(((\'A\',A0,sA,AFl,AFg,0),(\'B\',B0,sB,BLl,BLg,1))):\n  m,mgg=mb.copy(),mg.copy()\n  for i in range(4):\n    if i>0:n,p=steps2[i-1];pm=mgg.copy();m=m+T(n,p,th,k);mgg=mgg+T(n,p,gt,k)\n    if i==0:a=fig.add_subplot(gs[3*blk,0]);b=fig.add_subplot(gs[3*blk+1,0],sharex=a);AX[blk]=(a,b)\n    else:a=fig.add_subplot(gs[3*blk:3*blk+2,i],sharey=AX[blk][0]);b=None\n    for x in ((a,b) if b else (a,)):x.axvspan(15,80,color=\'gold\',alpha=.18,zorder=0);x.axvspan(260,335,color=\'lightsteelblue\',alpha=.35,zorder=0)\n    a.errorbar(l,Y,yerr=sY,fmt=\'o\',color=\'k\',ms=4,alpha=.5,lw=1)\n    if i>0:a.plot(g,pm,color=\'grey\',ls=\'--\',lw=2,alpha=.6)\n    a.plot(g,mgg,color=cols[i],lw=3.5)\n    r=Y-m;c2=np.sum(r**2/sY**2);r2=np.corrcoef(Y,m)[0,1]**2;rms=np.sqrt(np.mean(r**2));stats[(nm,i)]=(rms,c2,r2)\n    a.text(.03,.04,f\'χ²={c2:.2f}   r²={r2:.3f}\',transform=a.transAxes,fontsize=19,color=\'0.15\',bbox=dict(fc=\'white\',ec=\'none\',alpha=.8))\n    if b is None:pass\n    else:b.errorbar(l,r,yerr=sY,fmt=\'o\',color=\'k\',ms=4,alpha=.5,lw=1);b.axhline(0,color=\'grey\',lw=1.2)\n    if b is not None and i<3:\n      n2,p2=steps2[i];b.plot(g,p2[2]*np.cos(n2*gt)+p2[3]*np.sin(n2*gt),color=cols[i+1],lw=3.5)\n      b.text(.03,.05,f\'RMS={rms:.3f}; next: n={n2}, C={np.hypot(p2[2],p2[3]):.3f}\',transform=b.transAxes,fontsize=17,color=\'0.15\',bbox=dict(fc=\'white\',ec=\'none\',alpha=.8))\n    elif b is not None:b.text(.03,.05,f\'RMS={rms:.3f}; final residual\',transform=b.transAxes,fontsize=17,color=\'0.15\',bbox=dict(fc=\'white\',ec=\'none\',alpha=.8))\n    if True:a.set_title(\'CoE model (before)\' if i==0 else (\'+ shared n=1\' if i==1 else (\'+ shared n=1–2\' if i==2 else \'after all: n=1–3\')),fontweight=\'bold\',fontsize=21)\n    if i==0:a.text(-.32,.5,nm,transform=a.transAxes,fontsize=34,fontweight=\'bold\',va=\'center\')\n    a.set_xlim(0,360);a.set_xticks([0,90,180,270,360])\n    if i>0:a.tick_params(labelleft=False);a.set_xticks([90,180,270,360])\n    if blk==0 and i>0:a.tick_params(labelbottom=False)\n    if blk==0 and b is not None:b.tick_params(labelbottom=False)\n    if i==0:a.tick_params(labelbottom=False)\n    if blk==1 or i==0:(b if b is not None else a).set_xlabel(\'ℓ [deg]\',fontsize=21) if blk==1 else None\n  AX[blk][0].set_ylabel(\'A/A$_{all}$\' if nm==\'A\' else \'B/B$_{all}$\',fontsize=22);AX[blk][1].set_ylabel(f\'{nm} residual\',fontsize=22)\nh=[Line2D([],[],color=\'k\',marker=\'o\',ls=\'\',alpha=.5,label=\'reconstruction, 313 clusters (±σ Migkas Fig. 6)\'),Line2D([],[],color=\'purple\',lw=3,label=\'CoE model (A: AF_smooth, B: sheet BL)\')]+[Line2D([],[],color=cols[j],lw=3,label=f\'shared n={j} (model / component)\') for j in (1,2,3)]+[Line2D([],[],color=\'grey\',ls=\'--\',lw=2,label=\'previous model\'),plt.Rectangle((0,0),1,1,color=\'gold\',alpha=.3,label=\'ℓ 15–80°\'),plt.Rectangle((0,0),1,1,color=\'lightsteelblue\',alpha=.5,label=\'ℓ 260–335°\')]\nfig.legend(handles=h,loc=\'lower center\',ncol=3,frameon=False,fontsize=19)\nfig.subplots_adjust(left=.07,right=.99,top=.96,bottom=.14,wspace=.06);fig.savefig(\'recon_n3_2x4_r1.png\',dpi=90,facecolor=\'white\')\nfor kk,v in stats.items():print(kk,\'RMS=%.3f chi2=%.2f r2=%.3f\'%v)\nfor n,p in steps2:print(\'n\',n,\'C=%.4f\'%np.hypot(p[2],p[3]),\'phase=%.1f\'%(np.degrees(np.arctan2(p[3],p[2]))%360),\'offA=%.4f offB=%.4f\'%(p[0],p[1]))\n'
exec(compile(SOURCE, '/workspace/cosmology/figures/fourier_resid/recon_n3_2x4_r1.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/cosmology/figures/fourier_resid/recon_n3_2x4_r1.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
