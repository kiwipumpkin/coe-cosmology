# PROVENANCE COPY of /workspace/redweb2/stochfield/deg.py (original generating script, rev6 figure; run in its original directory).
# Output figure copied here as L23_degeneracy.png. Not re-run for rev7: see Appendix B.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,sys
sys.argv=['fit.py','30','bg'];src=open('fit.py').read().split('rng=np.random')[0];exec(src)
C=json.load(open('chosen_l30_bg.json'));p0=C['p']
def c2(**kw):q=dict(p0,**kw);Df=F.Dl(q);return chi2(Df+Dfix,Df)[0]
# slice 1: thin shell (width 0.3 Gpc): shell centre vs amplitude
cw=np.linspace(0.2,4.4,36);lA=np.linspace(-11,-6,41);S1=np.array([[c2(chiw=c,sw=0.3,lgA=a) for c in cw] for a in lA])
# slice 2: slope vs amplitude at chosen window
ns=np.linspace(-2,2,33);lA2=np.linspace(-11,-7,41);S2=np.array([[c2(n=n,lgA=a) for n in ns] for a in lA2])
np.savez('deg_l30_bg.npz',cw=cw,lA=lA,S1=S1,ns=ns,lA2=lA2,S2=S2)
fig,(a,b)=plt.subplots(1,2,figsize=(13,5));m1=S1.min();m2=S2.min()
for ax,X,Y,S,m,xl in ((a,cw,lA,S1,m1,'shell centre χ_w (Gpc), width 0.3 Gpc'),(b,ns,lA2,S2,m2,'P_δ slope n')):
    cs=ax.contourf(X,Y,np.log10(S),30,cmap='viridis');plt.colorbar(cs,ax=ax,label='log10 χ² (Planck+CV)')
    ax.contour(X,Y,S-m,[2.3,6.18,11.8],colors=['w','w','w'],linestyles=['-','--',':']);i,j=np.unravel_index(np.argmin(S),S.shape);ax.plot(X[j],Y[i],'r*',ms=12,label=f'min χ²={m:.1f}')
    ax.set_xlabel(xl);ax.set_ylabel('log10 A (P_δ amplitude)');ax.legend(loc='lower center')
a.set_title('(d) radial–amplitude degeneracy (Δχ² 2.3/6.2/11.8 contours)');b.set_title('(d) slope–amplitude at chosen window')
fig.tight_layout();fig.savefig('d_degeneracy_l30_bg.png',dpi=100)
print('slice1 min',m1,'at',cw[np.argmin(S1.min(0))],'; per-centre best lgA:',[(round(c,2),round(lA[np.argmin(S1[:,j])],2),round(S1[:,j].min(),1)) for j,c in enumerate(cw)][::5])
