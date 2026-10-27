# AUTO-WRAP generator for mc613_ladder.png
"""mc613_ladder.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: realweb/mc/mc_samples.py (12 random MC webs, rng seed 99)
Embedded source: /workspace/coepdf/run_out/run9/realweb/mc/mc_samples.py  (verbatim below, unmodified)
Run:  python3 mc613_ladder.py   -> writes mc613_ladder_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/run_out/run9/realweb/mc so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'mc613_ladder' + '_regen.png')
SRC_CWD = '/workspace/coepdf/run_out/run9/realweb/mc'
MATCH = ['mc613_ladder.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_mc613_ladder_')
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
sys.argv = ['/workspace/coepdf/run_out/run9/realweb/mc/mc_samples.py'] + []

SOURCE = 'import numpy as np,json,matplotlib;matplotlib.use(\'Agg\');import matplotlib.pyplot as plt\nexec(open(\'/workspace/redweb2/realweb/analyze.py\').read().split(\'O=\')[0])\nS=np.load(\'mc_sky.npy\',mmap_mode=\'r\');sm=np.load(\'/workspace/redweb2/sm_red.npy\');RS=[0.5,0.9,2.0,4.0,8.0,14.0]\nt=np.loadtxt(\'/workspace/joint/lcdm.txt\');l=np.arange(2,31);Lt=np.interp(l,t[:,0],t[:,1])\nE=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];G=np.logspace(-3.5,0.3,60)\nks=np.random.default_rng(99).choice(len(S),12,replace=False);out={}\nfor k in ks:\n  k=int(k);fig,ax=plt.subplots(6,1,figsize=(11,30));out[k]={}\n  for a,iR,R in zip(ax,range(6),RS):\n    sk=np.array(S[k][iR]);cs=[chi(Dl(sm+e*sk)) for e in G];eb=G[int(np.argmin(cs))];cb=min(cs)\n    Db=Dl(sm+eb*sk)[2:31];y=Db[2:];slope=np.polyfit(l[2:],np.log(y),1)[0];over=Db[:3].max()/Pl[:3,1].max()\n    ok=cb<=100 and abs(slope)<0.05 and over<2.5\n    out[k][R]=dict(best_e=float(eb),chi2=float(cb),slope=float(slope),lowl_over=float(over),guess=\'PASS\' if ok else \'FAIL\')\n    a.plot(l,Lt,\'k-\',lw=2,label=\'ΛCDM χ²=20.1\');D0=Dl(sm);a.plot(l,D0[2:31],\'r--\',lw=2,label=f\'smooth red χ²={chi(D0):.1f}\')\n    for j,e in enumerate(E):\n      D=Dl(sm+e*sk);a.plot(l,D[2:31],color=plt.cm.viridis(j/9),label=f\'{e*100:g}% χ²={chi(D):.4g}\')\n    a.plot(l,Db,\'m-\',lw=2.5,label=f\'best {eb*100:.3g}% χ²={cb:.0f}\')\n    a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt=\'o\',color=\'k\',ms=4,label=\'Planck\')\n    a.set_yscale(\'log\');a.set_ylim(1e-3,None);a.set_ylabel(\'D_ℓ (μK²)\')\n    a.set_title(f\'red + MC web #{k}, kernel at {R} Gpc — my guess: {out[k][R]["guess"]}\',color=\'g\' if ok else \'r\')\n    a.legend(loc=\'lower center\',ncol=4,fontsize=7,framealpha=.85)\n  ax[-1].set_xlabel(\'ℓ\');plt.tight_layout();plt.savefig(f\'samples/mc{k}_ladder.png\',dpi=80);plt.close()\njson.dump(out,open(\'samples/guesses.json\',\'w\'),indent=1)\nfor k in out:print(k,{R:(v[\'guess\'],round(v[\'chi2\']),f"{v[\'best_e\']*100:.3g}%") for R,v in out[k].items()})\n'
exec(compile(SOURCE, '/workspace/coepdf/run_out/run9/realweb/mc/mc_samples.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/run_out/run9/realweb/mc/mc_samples.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
