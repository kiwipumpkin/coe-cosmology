# Level-set ladder for the acoustic free_W best fit: R=2,4,8,14 Gpc, window centre R-0.056, width 0.05; amplitude x0.01..x100.
# Field only (no deterministic background); exact l<=50, Limber l>50 (no exact-node correction -> l<150 overestimated ~1.5-5x for edge shells).
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.signal import find_peaks,savgol_filter
import acoustic_model as M
from acoustic_model import AFwd,Gc,rho,T0,LMAX
import acoustic_fit as A
CK='acoustic_ladder.json';S=json.load(open(CK)) if os.path.exists(CK) else {}
J=json.load(open('acoustic_fit.json'));p0=dict(J['free_W']['best']['p']);MULT=np.geomspace(0.01,100,9);ll=np.arange(LMAX+1)
def lim(F,p,l):
    a=max(0.01,p['chiw']-6*p['sw']);b=min(F.R,p['chiw']+6*p['sw']);ch=np.linspace(a,b,600);KM=np.interp(ch,F.ch,F.KM)
    W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);k=(l[:,None]+0.5)/ch[None];Pp=(4*np.pi*Gc*rho)**2*F.Pd(p,k)/k**4
    return np.trapezoid((KM*W)**2*k**2/3*Pp/ch**2,ch,axis=1)*l*(l+1)/(2*np.pi)*T0**2
def pk(D):s=savgol_filter(D[100:],61,3);q,_=find_peaks(s,prominence=0.02*s.max());return (q+100)[:5].tolist()
PlP=J['peaks_planck_binned_spline'];pl1=float(A.Pb[(A.Pb[:,0]>150)&(A.Pb[:,0]<300),1].max())
if __name__=="__main__":
 for R in (2.0,4.0,8.0,14.0):
     key=f'R{R:g}'
     if key in S:continue
     F=AFwd(R=R,lmax=M.LX,kmax=200,nk=500);p=dict(p0,chiw=R-0.056,sw=0.05);t=F.terms(p)
     D1=np.r_[0,0,t['scalar']+t['mlens'],lim(F,p,np.arange(M.LX+1,LMAX+1))]
     h=float(D1[180:300].max());lev=[float(m*h) for m in MULT];br=[(float(MULT[i]),float(MULT[i+1])) for i in range(8) if lev[i]<=pl1<=lev[i+1]]
     S[key]=dict(D1=D1.tolist(),peaks=pk(D1),first_peak_height_x1=h,bracket=br,mult_to_match=pl1/h);json.dump(S,open(CK,'w'));print(key,S[key]['peaks'],'x1 peak %.4g'%h,'bracket',br,'mult needed %.3g'%(pl1/h),flush=True)
 fig,axs=plt.subplots(2,2,figsize=(15,11));rng=np.random.default_rng(42);cs=plt.cm.viridis(np.linspace(0,1,9))
 for ax,R in zip(axs.ravel(),(2,4,8,14)):
     s=S[f'R{R:g}'];D1=np.array(s['D1']);l=ll[2:]
     ax.plot(l,A.Lt[2:],'k-',lw=1.8,label='ΛCDM');ax.plot(A.Pf[:,0],A.Pf[:,1],'.',color='grey',ms=1.5,alpha=.6,label='Planck 2018')
     for j,m in enumerate(MULT):
         Dm=m*D1[2:];Dr=Dm*np.random.default_rng(42+j).chisquare(2*l+1)/(2*l+1)
         ax.plot(l,Dm,'-',color=cs[j],lw=.6,alpha=.3);ax.plot(l,Dr,'-',color=cs[j],lw=.6,label=f'×{m:.3g}')
     ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(1,1e6);ax.set_xlim(2,2500);ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)')
     ax.set_title(f'R={R} Gpc, shell {R-0.056:.3f}±0.05 Gpc; peaks {s["peaks"]}; Planck-peak bracket {s["bracket"]}',fontsize=8);ax.legend(loc='lower center',ncol=6,fontsize=6,framealpha=.85)
 fig.suptitle(f'acoustic free_W ladder (one realization per level, seed 42+j; faint = mean). Planck peaks {PlP}',fontsize=10);fig.tight_layout();fig.savefig('acoustic_ladder.png',dpi=100)
