# AUTO-WRAP generator for acoustic_cl_fixed_shell_real.png
"""acoustic_cl_fixed_shell_real.png -- generating script (placed beside the PNG per the every-PNG-has-its-.py rule).

Provenance: byte-identical to run9/stochfield/acoustic_cl_fixed_shell_real.png, written by stochfield/real_figs.py
Embedded source: /workspace/coepdf/run_out/run9/stochfield/real_figs.py  (verbatim below, unmodified)
Run:  python3 acoustic_cl_fixed_shell_real.py   -> writes acoustic_cl_fixed_shell_real_regen.png beside this file (the original PNG is never overwritten).
The embedded code runs with cwd = /workspace/coepdf/run_out/run9/stochfield so its input files resolve; any other files it
tries to write are diverted to a scratch dir.
"""
import os, sys, builtins, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'acoustic_cl_fixed_shell_real' + '_regen.png')
SRC_CWD = '/workspace/coepdf/run_out/run9/stochfield'
MATCH = ['acoustic_cl_fixed_shell_real.png']       # basename(s) the source saves this figure under
SCRATCH = tempfile.mkdtemp(prefix='regen_acoustic_cl_fixed_shell_real_')
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
sys.argv = ['/workspace/coepdf/run_out/run9/stochfield/real_figs.py'] + []

SOURCE = '# One random sky per model (fixed seed 42): D_real = D_bg + D_field * chi2_{2l+1}/(2l+1); background is deterministic (no scatter).\nimport json,numpy as np,matplotlib;matplotlib.use(\'Agg\');import matplotlib.pyplot as plt\nimport acoustic_fit as A\nfrom acoustic_model import Dfield,Dfix,LMAX\nfrom fwd import Fwd\ndef real(Df,l,seed=42):return Df*np.random.default_rng(seed).chisquare(2*l+1)/(2*l+1)\nJ=json.load(open(\'acoustic_fit.json\'));Lt=A.Lt;Pf=A.Pf;Pb=A.Pb;ll=np.arange(LMAX+1)\nfor nm in (\'free_W\',\'fixed_shell\'):\n    S=J[nm];E=S[\'best\'];p=E[\'p\'];Df=Dfield(p)*np.array(S[\'corr\']);Dm=Df+Dfix;Dr=Dfix+real(Df,np.maximum(ll,1))\n    st={k:A.stats(Dr,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500)).items()};S[\'realization_seed42\']=st\n    fig,ax=plt.subplots(2,1,figsize=(12,10))\n    for a,(l0,l1),key in zip(ax,((2,30),(2,LMAX)),(\'l2_30\',\'l2_2500\')):\n        m=(ll>=l0)&(ll<=l1);a.plot(ll[m],Lt[m],\'k-\',lw=1.8,label=f\'ΛCDM  χ²(Planck err)={J["lcdm"][key]["chi2_planck"]:.1f}, Pearson r²={J["lcdm"][key]["r2"]:.3f}\')\n        cv=np.sqrt(2/(2*np.maximum(ll,1)+1))*Df;a.fill_between(ll[m],(Dm-cv)[m],(Dm+cv)[m],color=\'r\',alpha=.2,label=\'±1σ cosmic variance (field part)\')\n        a.plot(ll[m],Dm[m],\'r-\',lw=.7,alpha=.35,label=\'CoE ensemble mean (faint)\')\n        a.plot(ll[m],Dr[m],\'r-\',lw=1.0 if l1>30 else 1.6,marker=\'\' if l1>30 else \'o\',ms=3,label=f\'CoE one realization (seed 42, {nm})  χ²(Planck err)={st[key]["chi2_planck"]:.1f}, Pearson r²={st[key]["r2"]:.3f}\')\n        mm=(Pf[:,0]>=l0)&(Pf[:,0]<=l1)\n        if l1==30:a.errorbar(Pf[mm,0],Pf[mm,1],yerr=[Pf[mm,2],Pf[mm,3]],fmt=\'o\',color=\'grey\',ms=4,label=\'Planck 2018\');a.set_ylim(0,2500)\n        else:a.plot(Pf[mm,0],Pf[mm,1],\'.\',color=\'grey\',ms=1.5,alpha=.5,label=\'Planck 2018 (unbinned)\');a.set_xscale(\'log\');a.set_ylim(0,7000)\n        a.set_xlabel(\'ℓ\');a.set_ylabel(\'D_ℓ (μK²)\');a.legend(loc=\'lower center\' if l1==30 else \'upper left\',fontsize=7,framealpha=.85)\n    ax[0].set_title(f\'acoustic-feature stochastic field, R=4 Gpc, {nm}: one realization\');fig.tight_layout();fig.savefig(f\'acoustic_cl_{nm}_real.png\',dpi=100);plt.close()\n    print(nm,st)\njson.dump(J,open(\'acoustic_fit.json\',\'w\'),indent=1)\n# stage A (l<=30 / l<=100) panels with a realization\nfor tag in (\'l30_bg\',\'l30_nobg\',\'l100_bg\',\'l100_nobg\'):\n    R=json.load(open(f\'fit_{tag}.json\'));C=json.load(open(f\'chosen_{tag}.json\'));P_=np.array(R[\'planck\']);l=P_[:,0].astype(int);F=Fwd(lmax=int(R[\'LMAX\']))\n    Df=F.Dl(C[\'p\']);Dfx=np.array(R[\'Dfix\']);Dm=Df+Dfx;Dr=Dfx+real(Df,l)\n    e=np.where(Dr>P_[:,1],P_[:,3],P_[:,2]);c2=float(np.sum((Dr-P_[:,1])**2/e**2));r2=float(np.corrcoef(Dr,P_[:,1])[0,1]**2)\n    fig,ax=plt.subplots(figsize=(11,5));ax.plot(l,R[\'lcdmD\'],\'k-\',lw=2,label=f\'ΛCDM χ²(Planck err)={R["lcdm"]["chi2_planck"]:.1f}, Pearson r²={R["lcdm"]["r2"]:.3f}\')\n    cv=np.sqrt(2/(2*l+1))*Df;ax.fill_between(l,Dm-cv,Dm+cv,color=\'r\',alpha=.2,label=\'±1σ cosmic variance (field part)\')\n    ax.plot(l,Dm,\'C0-\',lw=.7,alpha=.35,label=\'ensemble mean (faint)\');ax.plot(l,Dr,\'C0-o\',ms=3,lw=1.3,label=f\'one realization (seed 42): χ²(Planck err)={c2:.1f}, Pearson r²={r2:.3f}\')\n    ax.errorbar(l,P_[:,1],yerr=[P_[:,2],P_[:,3]],fmt=\'o\',color=\'grey\',ms=3,label=\'Planck 2018\');ax.set_ylim(0,2500 if R[\'LMAX\']<=30 else 3300)\n    ax.set_xlabel(\'ℓ\');ax.set_ylabel(\'D_ℓ (μK²)\');ax.set_title(f\'stochastic field, R=4 Gpc, one realization [{tag}]\');ax.legend(loc=\'lower center\',fontsize=7,framealpha=.85);fig.tight_layout();fig.savefig(f\'a_cl_{tag}_real.png\',dpi=100);plt.close()\n    print(tag,\'real chi2\',round(c2,1),\'r2\',round(r2,3))\n'
exec(compile(SOURCE, '/workspace/coepdf/run_out/run9/stochfield/real_figs.py', 'exec'), {'__name__': '__main__', '__file__': '/workspace/coepdf/run_out/run9/stochfield/real_figs.py'})
print('wrote', OUT if os.path.exists(OUT) else '(no matching figure produced)')
