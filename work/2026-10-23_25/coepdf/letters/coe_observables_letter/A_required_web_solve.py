# Inverse problem: find sky perturbation delta(l,b) (+ = matter/clusters, - = voids) such that
# A_new = y + strip_avg(delta) (Migkas strips: |dl|<45 deg, all b, cos b weight, at the 71 bin centres),
# joined connect-the-dots (np.interp periodic), FFT, dipole removed, has power share per m equal to the
# Planck 2018 TT (2l+1)C_l share for m=l=2..35. Total m=2..35 power of y kept fixed (alt totals also run).
# Step 1: latitude-independent delta(l) (minimal norm for a strip average) on a 1 deg grid, optimised with a small norm penalty.
# Step 2: same Delta reproduced by a finite set of Gaussian clusters (+, 3-8 deg) and voids (-, 15-30 deg) via NNLS.
# Writes A_required_web_results.npz used by the A_required_web_*.py plot scripts.
import numpy as np
from scipy.optimize import least_squares,nnls
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
o=np.argsort(d['l_deg']);l=d['l_deg'][o];y=d['A_over_Aall'][o];elo=d['A_err_low'][o];ehi=d['A_err_high'][o]
N=2**14;g=np.arange(N)*360/N
def spec(v):return 4*np.abs(np.fft.rfft(np.interp(g,l,v,period=360))/N)**2
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0].astype(int);Pl=(2*L+1)*2*np.pi*p[:,1]/(L*(L+1))
M=35;m=np.arange(2,M+1);pl=Pl[(L>=2)&(L<=M)];pl=pl/pl.sum()
# strip operator on latitude-independent delta(l), 1 deg grid
lg=np.arange(0.5,360,1.0)
S=np.array([(np.abs((lg-c+180)%360-180)<45).astype(float) for c in l]);S/=S.sum(1,keepdims=True)
P0=spec(y);T0=P0[2:M+1].sum()
def run(Tscale,lam=1e-3):
  T=T0*Tscale
  def res(x):
    D=S@x;P=spec(y+D)[2:M+1];return np.concatenate([(P-T*pl)/T, lam*x/np.sqrt(len(x))])
  # start: y's phases, Planck amplitudes, sampled, divided by nothing (optimiser handles strip)
  r=least_squares(res,np.zeros(360),method='trf',max_nfev=4000,x_scale=0.05)
  return r.x
out={}
for ts in [1.0,0.5,2.0]:
  x=run(ts);D=S@x;P=spec(y+D);a=P[2:M+1]/P[2:M+1].sum()
  out[ts]=(x,D,a);print('Tscale',ts,'rms resid after',np.sqrt(((a-pl)**2).mean()),'delta rms',x.std(),'peak',np.abs(x).max(),'Delta rms',D.std())
x,D,a=out[1.0]
a0=P0[2:M+1]/T0;print('rms resid before',np.sqrt(((a0-pl)**2).mean()))
# transfer function of 90 deg boxcar
H=lambda mm:np.sinc(mm*45/180)
# --- finite clusters/voids via NNLS on the 2 deg sphere ---
lg2=np.radians(np.arange(1,360,2));bg2=np.radians(np.arange(-89,90,2));LL,BB=np.meshgrid(lg2,bg2);W=np.cos(BB)
xyz=np.stack([np.cos(BB)*np.cos(LL),np.cos(BB)*np.sin(LL),np.sin(BB)]).reshape(3,-1)
St=np.array([((np.abs(np.angle(np.exp(1j*(LL-np.radians(c)))))<np.radians(45))*W).ravel() for c in l]);St/=St.sum(1,keepdims=True)
cols=[];meta=[]
for kind,sign,sigs in [('cluster',1,[3,5,8]),('void',-1,[15,22,30])]:
  for sg in sigs:
    for lc in np.arange(0,360,10):
      for bc in [-60,-30,0,30,60]:
        v=np.array([np.cos(np.radians(bc))*np.cos(np.radians(lc)),np.cos(np.radians(bc))*np.sin(np.radians(lc)),np.sin(np.radians(bc))])
        bl=np.exp(-0.5*(np.arccos(np.clip(v@xyz,-1,1))/np.radians(sg))**2)
        cols.append(sign*(St@bl));meta.append((kind,sg,lc,bc,sign))
C=np.array(cols).T
# allow free constant (overall normalisation) via +/- columns
C2=np.hstack([C,np.ones((71,1)),-np.ones((71,1))])
mu=3e-4  # small ridge so NNLS avoids huge cancelling blobs
w,rn=nnls(np.vstack([C2,mu*np.eye(C2.shape[1])]),np.concatenate([D,np.zeros(C2.shape[1])]),maxiter=50000)
Dfit=C2@w;meta=np.array(meta,dtype=object);wb=w[:-2]
keep=wb>1e-3*wb.max()
Pf=spec(y+Dfit);af=Pf[2:M+1]/Pf[2:M+1].sum()
print('NNLS fit rms of Delta mismatch',np.sqrt(((Dfit-D)**2).mean()),'resid after (blobs)',np.sqrt(((af-pl)**2).mean()),'const',w[-2]-w[-1])
for kind in ['cluster','void']:
  k=keep&(meta[:,0]==kind);print(kind,'n',k.sum(),'peak amps min/median/max',wb[k].min() if k.any() else 0,np.median(wb[k]) if k.any() else 0,wb[k].max() if k.any() else 0)
mapb=(np.array([ (np.exp(-0.5*(np.arccos(np.clip(np.array([np.cos(np.radians(b))*np.cos(np.radians(lc)),np.cos(np.radians(b))*np.sin(np.radians(lc)),np.sin(np.radians(b))])@xyz,-1,1))/np.radians(sg))**2)) for (_,sg,lc,b,_) in meta[keep]]).T*(wb[keep]*meta[keep][:,4].astype(float))).sum(1).reshape(LL.shape)
np.savez('A_required_web_results.npz',l=l,y=y,elo=elo,ehi=ehi,lg=lg,x=x,D=D,a=a,a0=a0,pl=pl,m=m,
 x05=out[0.5][0],a05=out[0.5][2],x2=out[2.0][0],a2=out[2.0][2],Dfit=Dfit,af=af,mapb=mapb,
 keep_meta=np.array([list(r)[:4] for r in meta[keep]],dtype=object),keep_w=wb[keep],P0=P0,Pnew=spec(y+D),Pblob=Pf,Pl=Pl,L=L,allow_pickle=True)
# m multiples of 4: share achievable?
for mm in [4,8,12,16]:print('m',mm,'before',a0[mm-2],'after',a[mm-2],'planck',pl[mm-2])
