# AUTO-WRAP generator for L11_web_slab.png
"""L11_web_slab.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to cosmology-course/.../constrained/red_web_s2000_bodies_1pct.png; generator recovered as coepdf/book/recovered/f06.py
Embedded source: /workspace/coepdf/book/recovered/f06.py  (verbatim below, unmodified)
Run:  python3 L11_web_slab.py   -> writes L11_web_slab_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/book/recovered so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'L11_web_slab' + '_regen.png')
SRC_CWD = '/workspace/coepdf/book/recovered'
MATCH = ['red_web_s2000_bodies_1pct.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_L11_web_slab_')
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
sys.argv = ['/workspace/coepdf/book/recovered/f06.py'] + []

SOURCE = '"""Plain side view of one generated network: slab |z| < half-thickness, void core spheres cut by the slab,\nclusters, and every void-cluster link with both ends in the slab. Full region + zoom."""\nimport os,sys,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt\nfrom matplotlib.collections import LineCollection,PatchCollection\nME=\'/workspace/coe-cosmology-gh/draft/current/analyses/void_mass_network\'; sys.path.insert(0,ME); import vmn_gen as G\nMmin=float(sys.argv[1]) if len(sys.argv)>1 else 3e15\ng=G.generate(Mmin=Mmin,seed=2000,comp_mode="lcdm"); S=g["S"]/1000; XC=g["XC"]/1000; Rv=g["Rv"]/1000; L=g["links"]; R=g["params"]["R_out"]/1000\ndef draw(ax,lim,hz,cen=(0,0)):\n    sv=(np.abs(S[:,2])<hz)&(np.abs(S[:,0]-cen[0])<lim+2)&(np.abs(S[:,1]-cen[1])<lim+2)\n    rc=np.sqrt(np.maximum(Rv[sv]**2-np.maximum(np.abs(S[sv,2])-0,0)**2,0))      # radius of the slice through z=0\n    pc=PatchCollection([plt.Circle(S[i,:2],r) for i,r in zip(np.where(sv)[0],rc) if r>0],facecolor="#2471a3",alpha=0.35,edgecolor="none",zorder=1)\n    ax.add_collection(pc)\n    sc=(np.abs(XC[:,2])<hz)&(np.abs(XC[:,0]-cen[0])<lim+2)&(np.abs(XC[:,1]-cen[1])<lim+2)\n    lk=L[sc[L[:,0]]&sv[L[:,1]]]\n    seg=np.stack([XC[lk[:,0],:2],S[lk[:,1],:2]],1)\n    ax.add_collection(LineCollection(seg,colors="0.3",linewidths=0.35 if lim>8 else 0.8,alpha=0.6,zorder=2))\n    m=0.01*g["MC"][sc]; ax.scatter(XC[sc,0],XC[sc,1],s=(2 if lim>8 else 8)+(30 if lim>8 else 80)*(m/m.max())**1.5,c="#c0392b",edgecolors="none",zorder=3); ax.text(0.02,0.02,f"cluster masses at 1%%: {m.min():.2e} to {m.max():.2e} Msun (dot size relative)",transform=ax.transAxes,fontsize=10,bbox=dict(fc="w"))\n    ax.add_patch(plt.Circle((0,0),14.3,fill=False,color="k",lw=1.4,zorder=4)); ax.plot(0,0,"k+",ms=12,mew=2,zorder=5)\n    ax.add_patch(plt.Circle((0,0),0.5,fill=False,color="k",lw=0.8,ls=":",zorder=4))\n    for rr,cc in ((0.23,"m"),(1,"tab:orange"),(2,"g"),(3,"c"),(4,"tab:purple")): ax.add_patch(plt.Circle((0,0),rr,fill=False,color=cc,lw=1.2,ls="--",zorder=4))\n    ax.set_xlim(cen[0]-lim,cen[0]+lim); ax.set_ylim(cen[1]-lim,cen[1]+lim); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")\n    return sv.sum(),sc.sum(),len(lk)\nfig,ax=plt.subplots(2,1,figsize=(10,20.5))\na,b,c=draw(ax[0],R+0.5,0.4); ax[0].set_title(f"Web used on red (seed 2000, lcdm) at 1% mass (full-mass cutoff {Mmin:.0e} Msun, i.e. {0.01*Mmin:.0e} at 1%): slab |z| < 0.4 Gpc of the {R:.0f} Gpc region\\n{b} clusters, {c} void–cluster links shown; black circle = our 14.3 Gpc sphere, + = us",fontsize=12)\na,b,c=draw(ax[1],2.5,0.4,cen=(0,0)); ax[1].set_title(f"Zoom (5 × 5 Gpc) around us; dashed = 0.23/1/2/3/4 Gpc evaluation radii: blue = void cores (non-overlapping),\\nred = clusters (size ~ mass), grey = links from each cluster to every void it borders",fontsize=12)\nfig.text(0.5,0.008,"Blue discs = core spheres of voids centred in the slab (each void is the whole cell between its clusters).\\nRed dots = clusters at web nodes. Grey lines = bipartite links (void ↔ cluster only),\\ndrawn when both ends are in the slab. Dotted circle = 0.5 Gpc zone kept empty around us.",ha="center",fontsize=10.5)\nfig.tight_layout(rect=(0,0.045,1,1)); fig.savefig("/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained/red_web_s2000_bodies_1pct.png",dpi=110)\n'
exec(compile(SOURCE, '/workspace/coepdf/book/recovered/f06.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/book/recovered/f06.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
