# AUTO-WRAP generator for sky_mapB_al_remade.png
"""sky_mapB_al_remade.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: per sky.SOURCE.txt: /tmp r45c.py then manual crop to map B axes (crop step not scripted)
Embedded source: /tmp/scr/l50/r45c.py  (verbatim below, unmodified)
Run:  python3 sky_mapB_al_remade.py   -> writes sky_mapB_al_remade_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /tmp/scr/l50 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'sky_mapB_al_remade' + '_regen.png')
SRC_CWD = '/tmp/scr/l50'
MATCH = ['overlay_lift15_r45_bluered.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_sky_mapB_al_remade_')
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
sys.argv = ['/tmp/scr/l50/r45c.py'] + []

SOURCE = 'import json,numpy as np,matplotlib;matplotlib.use(\'Agg\');import matplotlib.pyplot as plt\nfrom matplotlib.lines import Line2D\nexec(open(\'/tmp/spots.py\').read());MM=json.load(open(\'/tmp/M25.json\'))\nZ={\'A\':[22,99,262,310],\'B\':[14,104,266,321]}\nTILT,ROLL=-15.0,45.0\ndef vec(l,b):l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)\ndef rot(ax,a):\n  ax=ax/np.linalg.norm(ax);a=np.radians(a);K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]])\n  return np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K\nR1=rot(np.array([0,1.,0]),TILT)           # CoE (180,0) -> b=+15\ncoe=R1@vec(180,0);R=rot(coe,ROLL)@R1      # +80 right-handed about CoE axis\ndef lm(l,b):m=vec(l,b)@R;return np.degrees(np.arctan2(m[...,1],m[...,0]))%360   # model-frame l\ndef moll(l,b):\n  lam=np.radians(((np.asarray(l,float)+180)%360)-180);b=np.radians(np.asarray(b,float));th=b.copy()\n  for _ in range(60):th=th-(2*th+np.sin(2*th)-np.pi*np.sin(b))/(2+2*np.cos(2*th)+1e-12)\n  return -2*np.sqrt(2)/np.pi*lam*np.cos(th),np.sqrt(2)*np.sin(th)\n# inverse mollweide grid\nX,Y=np.meshgrid(np.linspace(-2*np.sqrt(2),2*np.sqrt(2),1400),np.linspace(-np.sqrt(2),np.sqrt(2),700))\nok=(X/(2*np.sqrt(2)))**2+(Y/np.sqrt(2))**2<=1\nth=np.arcsin(np.clip(Y/np.sqrt(2),-1,1));B=np.degrees(np.arcsin(np.clip((2*th+np.sin(2*th))/np.pi,-1,1)))\nL=np.degrees(-np.pi*X/(2*np.sqrt(2)*np.cos(th)+1e-12))%360\nLM=lm(L,B);rows=[]\nfig,axs=plt.subplots(2,1,figsize=(15,15),facecolor=\'white\')\nfor ax,K in zip(axs,\'AB\'):\n  v=np.interp(LM,MM[\'l\'],MM[K],period=360);gm=np.abs(v[ok]).max()\n  img=np.ones(X.shape+(4,));img[...,3]=0\n  f=np.clip((v/gm+1)/2,0,1)\n  img[...,0]=f;img[...,1]=0;img[...,2]=1-f;img[...,3]=np.where(ok,.6,0)\n  ax.imshow(img,extent=(X.min(),X.max(),Y.min(),Y.max()),origin=\'lower\',interpolation=\'bilinear\')\n  for z in Z[K]:\n    d=(LM-z+180)%360-180;ax.contour(X,Y,np.where(ok,d,np.nan),[0],colors=\'k\',linewidths=2)\n  t=np.linspace(0,2*np.pi,400);ax.plot(2*np.sqrt(2)*np.cos(t),np.sqrt(2)*np.sin(t),\'k\',lw=1.5)\n  for lg in (0,90,270):ax.text(moll(lg,0)[0],-1.55,f\'{lg}°\',ha=\'center\',fontsize=12)\n  for nm_,(cl,cb),mk in [(\'CoE\',(180,0),\'X\')]:\n    c=R1@vec(cl,cb);c2=np.degrees([np.arctan2(c[1],c[0])%(2*np.pi),np.arcsin(c[2])]);xs,ys=moll(*c2);ax.plot(xs,ys,\'kX\',ms=16)\n  MK=[\'o\',\'s\',\'D\',\'^\',\'v\',\'p\',\'h\'];h=[Line2D([],[],ls=\'\',marker=\'X\',color=\'k\',ms=13,label=\'CoE (b = −15°; anti-CoE lifted +15°)\')];sc=0\n  for i,(nm,sl,sb,c) in enumerate(SPOTS[K]):\n    val=float(np.interp(lm(sl,sb),MM[\'l\'],MM[K],period=360));okk=(val>0)==(c==\'hot\');sc+=okk\n    rows.append((K,nm,sl,sb,round(float(lm(sl,sb)),1),round(val,2),okk))\n    xs,ys=moll(sl,sb);xs,ys=float(xs),float(ys)\n    ax.plot(xs,ys,MK[i],ms=28,mfc=CC[c],mec=\'k\',mew=2,zorder=7);ax.text(xs,ys,str(i+1),ha=\'center\',va=\'center\',color=\'white\',fontsize=15,fontweight=\'bold\',zorder=8)\n    ax.annotate(f"{nm} {val:+.1f}% {\'✓\' if okk else (\'?\' if abs(val)<0.1 else \'✗\')}",(xs,ys),xytext=(-20,14) if xs>0 else (20,14),textcoords=\'offset points\',ha=\'right\' if xs>0 else \'left\',fontsize=12.5,fontweight=\'bold\',color=CC[c],zorder=9,bbox=dict(fc=\'white\',ec=CC[c],boxstyle=\'round,pad=.25\'))\n    h.append(Line2D([],[],ls=\'\',marker=MK[i],mfc=CC[c],mec=\'k\',ms=13,label=f\'{i+1} {nm} ({sl}°, {sb:+g}°)\'))\n  h+=[plt.Rectangle((0,0),1,1,fc=\'red\',alpha=.45,label=f\'{K} residual high\'),plt.Rectangle((0,0),1,1,fc=\'blue\',alpha=.45,label=f\'{K} residual low\'),Line2D([],[],color=\'k\',lw=2,label=\'zero crossings (rotated)\')]\n  ax.legend(handles=h,loc=\'upper center\',bbox_to_anchor=(.5,-.02),ncol=3,fontsize=11.5,frameon=False)\n  ax.set_title(f\'{K} residual: axis lifted 15° , then 45° about CoE axis   —   {sc}/{len(SPOTS[K])} match\'+(\' (1 ?)\' if K==\'A\' else \'\'),fontsize=15)\n  ax.set_aspect(\'equal\');ax.axis(\'off\');ax.set_xlim(-2.9,2.9);ax.set_ylim(-1.65,1.5)\nfig.suptitle(\'Galactic Mollweide, ℓ = 0 at centre, ℓ increases to the left; spots at published (ℓ, b)\',fontsize=13,y=.995)\nfig.tight_layout();fig.savefig(\'overlay_lift15_r45_bluered.png\',dpi=110,facecolor=\'white\')\nfor r in rows:print(r)\n'
exec(compile(SOURCE, '/tmp/scr/l50/r45c.py', 'exec'), {'__name__': '__main__', '__file__': '/tmp/scr/l50/r45c.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
