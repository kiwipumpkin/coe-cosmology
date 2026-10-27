# AUTO-WRAP generator for top_new.png
"""top_new.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: /tmp/scr/l50/topnew.py (also coepdf/book/recovered/tmp_rescue/topnew.py); byte-identical copy /tmp/top_new.png
Embedded source: /tmp/scr/l50/topnew.py  (verbatim below, unmodified)
Run:  python3 top_new.py   -> writes top_new_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /tmp/scr/l50 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'top_new' + '_regen.png')
SRC_CWD = '/tmp/scr/l50'
MATCH = ['top_new.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_top_new_')
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
sys.argv = ['/tmp/scr/l50/topnew.py'] + []

SOURCE = 'src=open(\'/workspace/cosmology/figures/l50/boot_spots.py\').read().split("M=json.load")[0]\nexec(src)\nimport matplotlib;matplotlib.use(\'Agg\');import matplotlib.pyplot as plt\nfrom matplotlib.lines import Line2D\nexec(open(\'/tmp/spots.py\').read())\nZ={\'A\':[22,99,262,310],\'B\':[14,104,266,321]}\nfig,axs=plt.subplots(2,1,figsize=(18,11),facecolor=\'white\')\ncols={15:\'#3b0f4f\',25:\'#2c6e8f\',35:\'#3cbc63\'}\nfor ax,K in zip(axs,\'AB\'):\n  r25=base[K][25];gm=np.abs(r25).max();ll=np.arange(0,360.5,.5)\n  v=np.interp(ll,l,r25,period=360);f=np.clip((v/gm+1)/2,0,1)\n  img=np.stack([.3*f+.7,np.full_like(f,.7),.3*(1-f)+.7],-1)[None]\n  ax.imshow(img,extent=(0,360,-100,100),aspect=\'auto\',zorder=0)\n  h=[]\n  for w in (15,25,35):\n    y=base[K][w];o=np.argsort(l);ax.plot(np.r_[l[o],l[o][0]+360],np.r_[y[o],y[o][0]],color=cols[w],lw=2.6 if w==15 else 2,zorder=3)\n    imx,imn=np.argmax(y),np.argmin(y)\n    h.append(Line2D([],[],color=cols[w],lw=2.5,label=f\'±{w}°  (max {y[imx]:+.1f}% at {l[imx]:.0f}°, min {y[imn]:+.1f}% at {l[imn]:.0f}°)\'))\n  ax.axhline(0,color=\'grey\',lw=1)\n  for z in Z[K]:ax.axvline(z,color=\'k\',lw=1.5,ls=\'--\',zorder=2)\n  lo,hi=ax.get_ylim() if False else (min(base[K][15].min(),-1)*1.15,max(base[K][15].max(),1)*1.45)\n  ax.set_ylim(lo,hi)\n  MK=[\'o\',\'s\',\'D\',\'^\',\'v\',\'p\',\'h\']\n  for i,(nm,sl,sb,c) in enumerate(SPOTS[K]):\n    lm_={\'CMB dipole cold\':97.2,\'Bulk flow away\':113.9,\'CMB Cold Spot\':156.2,\'CMB dipole hot\':277.2,\'Bulk flow toward\':293.9,\'NVSS cold\':58.5,\'Hemi. asym. cold\':29.9,\'CatWISE cold\':72.2,\'NVSS hot\':238.5,\'Hemi. asym. hot\':209.9,\'CatWISE hot\':252.2,\'Loop I / NPS\':335.1}[nm]\n    ax.axvline(lm_,color=CC[c],lw=2,ls=\'-.\',zorder=2)\n    ax.plot(lm_,hi*.88,MK[i],ms=17,mfc=CC[c],mec=\'k\',mew=1.5,zorder=5);ax.text(lm_,hi*.88,str(i+1),ha=\'center\',va=\'center\',color=\'white\',fontsize=10,fontweight=\'bold\',zorder=6)\n    h.append(Line2D([],[],ls=\'\',marker=MK[i],mfc=CC[c],mec=\'k\',ms=11,label=f\'{i+1} {nm}\'))\n  h+=[plt.Rectangle((0,0),1,1,fc=(1,.7,.7),label=\'residual high\'),plt.Rectangle((0,0),1,1,fc=(.7,.7,1),label=\'residual low\'),Line2D([],[],color=\'k\',ls=\'--\',label=\'zero crossings\')]\n  ax.legend(handles=h,loc=\'upper center\',bbox_to_anchor=(.5,-.13),ncol=4,fontsize=11,frameon=False)\n  ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.set_xlabel(\'model-frame angle ℓ_m [deg]\',fontsize=13)\n  ax.set_ylabel(f\'{K} running mean residual\\n% of {K} (after n=1–3)\',fontsize=13)\nfig.tight_layout();fig.savefig(\'top_new.png\',dpi=110,facecolor=\'white\')\n'
exec(compile(SOURCE, '/tmp/scr/l50/topnew.py', 'exec'), {'__name__': '__main__', '__file__': '/tmp/scr/l50/topnew.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
