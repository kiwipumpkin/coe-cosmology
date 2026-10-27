# AUTO-WRAP generator for L13_target_profiles.png
"""L13_target_profiles.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: original f12 script never found; R13_target_profiles.py is the computed remake
Embedded source: /workspace/coepdf/fig7/R13_target_profiles.py  (verbatim below, unmodified)
Run:  python3 L13_target_profiles.py   -> writes L13_target_profiles_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/fig7 so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'L13_target_profiles' + '_regen.png')
SRC_CWD = '/workspace/coepdf/fig7'
MATCH = ['R13_target_profiles.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_L13_target_profiles_')
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
sys.argv = ['/workspace/coepdf/fig7/R13_target_profiles.py'] + []

SOURCE = '# rev7 remake of rev6 Fig 13 (target angular profiles matched by masses). Original script not found; rebuilt from the caption\'s method.\n# Dimensionless in R (lengths / R_ls), so it does not depend on D. Potential on the sphere |n|=1:\n#   Phi(theta) = M0/|n - 1.55 z| + sum_j m_j / |n - p_j|, M0 fixed (=1), p_j = 40 off-sphere sites (radii 0.5-0.9 and 1.15-1.4, polar angles 0-180 deg\n#   in the x-z plane, axisymmetrised by averaging over azimuth), m_j signed (clump>0, void<0) fitted by ridge least squares so that\n#   Phi(theta)-mean matches a*T(theta) for T = cos(theta/2), cos(1.5 theta), cos(2.5 theta), cos(2 theta).\nimport numpy as np,json,matplotlib;matplotlib.use(\'Agg\');import matplotlib.pyplot as plt\nth=np.radians(np.linspace(0,180,181));ph=np.linspace(0,2*np.pi,72,endpoint=False)\ndef ring(r,a):  # azimuth-averaged 1/|n-p| for p at radius r, polar angle a (an axisymmetric ring of mass)\n    n=np.stack([np.sin(th)[:,None]*np.cos(ph),np.sin(th)[:,None]*np.sin(ph),np.cos(th)[:,None]+0*ph],-1);p=np.array([r*np.sin(a),0,r*np.cos(a)])\n    p2=np.stack([r*np.sin(a)*np.cos(ph),r*np.sin(a)*np.sin(ph),r*np.cos(a)+0*ph],-1)\n    return np.mean(1/np.linalg.norm(n[:,None,:,:]-p2[None,:,None,:],axis=-1),axis=(1,2))\nsites=[(r,a) for r in (0.5,0.7,0.9,1.15,1.4) for a in np.radians(np.linspace(0,180,8))]\nG=np.stack([ring(r,a) for r,a in sites],1);G-=G.mean(0);P0=ring(1.55,0.0);P0-=P0.mean()\nT={\'half wave cos(θ/2)\':np.cos(th/2),\'1½ waves cos(1.5θ)\':np.cos(1.5*th),\'2½ waves cos(2.5θ)\':np.cos(2.5*th),\'two waves cos(2θ)\':np.cos(2*th)}\nfig,ax=plt.subplots(2,4,figsize=(22,9));out={}\nfor i,(k,t) in enumerate(T.items()):\n    t=t-t.mean();a=np.ptp(P0)/np.ptp(t);rhs=a*t-P0;lam=1e-6*np.trace(G.T@G)\n    m=np.linalg.solve(G.T@G+lam*np.eye(G.shape[1]),G.T@rhs);fit=P0+G@m;r2=1-np.sum((fit-a*t)**2)/np.sum((a*t)**2)\n    out[k]=dict(r2=float(r2),n_clump=int((m>0).sum()),n_void=int((m<0).sum()),sum_m=float(m.sum()))\n    A=ax[0,i];A.plot(np.degrees(th),a*t,\'k\',lw=2.5,label=\'target\');A.plot(np.degrees(th),fit,\'r--\',lw=2,label=f\'fit, R²={r2:.4f}\');A.plot(np.degrees(th),P0,\'0.6\',label=\'fixed mass alone\')\n    A.set_title(k);A.set_xlabel(\'θ from mass axis [deg]\');A.legend(fontsize=8)\n    B=ax[1,i];xs=np.array([r*np.sin(a_) for r,a_ in sites]);zs=np.array([r*np.cos(a_) for r,a_ in sites])\n    B.scatter(xs,zs,s=20+300*np.abs(m)/np.abs(m).max(),c=np.where(m>0,\'crimson\',\'royalblue\'),alpha=.7);B.plot(0,1.55,\'k*\',ms=14)\n    B.add_patch(plt.Circle((0,0),1,fill=False,ls=\'--\'));B.set_aspect(\'equal\');B.set_xlim(-.2,1.6);B.set_ylim(-1.6,1.7);B.set_xlabel(\'ρ / R\');B.set_ylabel(\'z / R\')\n    B.set_title(f\'{out[k]["n_clump"]} clump rings (red), {out[k]["n_void"]} void rings (blue); ★ fixed mass at 1.55R\',fontsize=8)\nfig.tight_layout();fig.savefig(__file__.replace(\'.py\',\'.png\'),dpi=90);json.dump(out,open(__file__.replace(\'.py\',\'.json\'),\'w\'),indent=1);print(out)\n'
exec(compile(SOURCE, '/workspace/coepdf/fig7/R13_target_profiles.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/fig7/R13_target_profiles.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
