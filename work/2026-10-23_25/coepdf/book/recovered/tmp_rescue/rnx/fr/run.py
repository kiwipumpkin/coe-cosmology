import numpy as np,csv,json,sys
sys.path.insert(0,'/workspace/cosmology/figures')
from make_A_vs_angle_axisymmetric_before_rotation import read_migkas,periodic_interp,CSV_AF
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import chi2 as C2, pearsonr
m=read_migkas(__import__('pathlib').Path('/workspace/cosmology/migkas/fig6_digitized.csv'))
rows=list(csv.DictReader(open(CSV_AF)))
x=np.array([float(r['ell_deg']) for r in rows]); af=np.array([float(r['AF_smooth']) for r in rows])
l=m['ell']; mod=periodic_interp(x,af,l); r=m['A']-mod
s=0.5*(m['err_low']+m['err_high']); w=1/s**2; N=len(l); print('N',N)
th=np.radians(l)
def design(n): return np.column_stack([np.ones(N)]+[f(k*th) for k in range(1,n+1) for f in (np.cos,np.sin)])
def fit(X):
    A=X*np.sqrt(w)[:,None]; b=r*np.sqrt(w); cov=np.linalg.inv(A.T@A); p=cov@(A.T@b)
    return p,cov,float(np.sum(w*(r-X@p)**2))
chi_raw=float(np.sum(w*r**2)); p0,_,chi0=fit(design(0))
out={'N':N,'chi2_raw':chi_raw,'offset':p0[0],'chi2_offset':chi0}; tab=[]; prev=chi0
for n in range(1,7):
    p,cov,c=fit(design(n)); tab.append(dict(n=n,chi2_cum=c,dchi2=prev-c,dof=N-(2*n+1),p_add=float(C2.sf(prev-c,2)))); prev=c
pF,covF,cF=fit(design(6))
for i,t in enumerate(tab):
    k=t['n']; a,b=pF[2*k-1],pF[2*k]; va,vb,cab=covF[2*k-1,2*k-1],covF[2*k,2*k],covF[2*k-1,2*k]
    C=np.hypot(a,b); sC=np.sqrt((a*a*va+b*b*vb+2*a*b*cab))/C
    ph=np.degrees(np.arctan2(b,a))%360; sph=np.degrees(np.sqrt(b*b*va+a*a*vb-2*a*b*cab)/C**2)/k
    # zero-crossing phase in sin form: C sin(k(l-l0))... report cos-phase phi_k: C cos(k l - phi_k)
    t.update(a=a,b=b,C=C,sC=sC,SNR=C/sC,phi_cos_deg=ph,sphi=sph,sin_phase_deg=((ph-90)%360)/k, peak_deg=ph/k)
out['table']=tab; out['chi2_full6']=cF
# comparison with rotation term actually used
Ar=np.array([float(q['A_rotated_smooth']) for q in rows]); rot=periodic_interp(x,Ar-af,l)
out['pearson_r2_resid_vs_rotationterm']=pearsonr(r,rot)[0]**2
out['chi2_after_rotation_term']=float(np.sum(w*(r-rot)**2))
p2,_,c2=fit(design(2)); rec2=design(2)@p2
out['pearson_r2_resid_vs_n1n2']=pearsonr(r,rec2)[0]**2
p1,_,c1=fit(design(1)); out['pearson_r2_resid_vs_n1']=pearsonr(r,design(1)@p1)[0]**2
json.dump(out,open('fourier_residuals.json','w'),indent=1,default=float)
with open('residuals.csv','w') as f:
    f.write('ell_deg,A_migkas,A_axisym_AF_smooth,residual,sigma\n')
    for v in zip(l,m['A'],mod,r,s): f.write(','.join(f'{q:.6g}' for q in v)+'\n')
plt.rcParams.update({'font.size':16,'figure.facecolor':'white'})
g=np.linspace(0,360,721); gt=np.radians(g)
G=lambda n:np.column_stack([np.ones_like(g)]+[f(k*gt) for k in range(1,n+1) for f in (np.cos,np.sin)])
fig,ax=plt.subplots(figsize=(13,7))
ax.errorbar(l,r,yerr=[m['err_low'],m['err_high']],fmt='o',color='k',ms=5,capsize=2,label='Migkas A/A_all − axisymmetric model (AF_smooth)')
ax.plot(g,G(1)@p1,color='tab:blue',lw=2.5,label=f'Fourier n≤1 (χ²={c1:.1f})')
ax.plot(g,G(2)@p2,color='tab:red',lw=2.5,label=f'Fourier n≤2 (χ²={c2:.1f})')
ax.plot(g,np.interp(g,x,Ar-af),'--',color='tab:green',lw=2,label='rotation term used in deck (Amp −0.06, φ=0.35 rad)')
ax.axhline(0,color='grey',lw=1); ax.set_xlim(0,360); ax.set_xlabel('Galactic longitude ℓ [deg]'); ax.set_ylabel('residual r(ℓ)')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.13),ncol=2,fontsize=13,frameon=False); fig.tight_layout(); fig.savefig('residuals_fourier_fit.png',dpi=200)
fig,ax=plt.subplots(1,2,figsize=(14,6))
n=[t['n'] for t in tab]
ax[0].bar(n,[t['C'] for t in tab],yerr=[t['sC'] for t in tab],color='tab:blue',capsize=5); ax[0].set_xlabel('harmonic n'); ax[0].set_ylabel('amplitude C_n'); ax[0].set_title('Amplitude per harmonic')
ax[1].bar(n,[t['dchi2'] for t in tab],color='tab:red'); ax[1].axhline(5.99,ls='--',color='k',label='Δχ²=5.99 (p=0.05, 2 dof)'); ax[1].set_xlabel('harmonic n'); ax[1].set_ylabel('Δχ² on adding n'); ax[1].set_title('Δχ² per harmonic')
ax[1].legend(loc='upper center',bbox_to_anchor=(0.5,-0.15),frameon=False); fig.tight_layout(); fig.savefig('harmonic_power.png',dpi=200)
print(json.dumps(out,indent=1,default=float))
