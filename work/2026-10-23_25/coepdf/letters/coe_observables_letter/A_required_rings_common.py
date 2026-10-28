# Helper (no PNG): rebuilds ring sky / strip average / spectra from A_required_rings_results.json, same definitions as A_required_rings_solve.py.
import numpy as np,json
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
o=np.argsort(d['l_deg']);l=d['l_deg'][o];y=d['A_over_Aall'][o];elo=d['A_err_low'][o];ehi=d['A_err_high'][o]
Ng=2**14;g=np.arange(Ng)*360/Ng
def spec(v):return 4*np.abs(np.fft.rfft(np.interp(g,l,v,period=360))/Ng)**2
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0].astype(int);Pl=(2*L+1)*2*np.pi*p[:,1]/(L*(L+1))
M=35;m=np.arange(2,M+1);pl=Pl[(L>=2)&(L<=M)];pl=pl/pl.sum()
lgd=np.arange(1.5,360,3);bgd=np.arange(-88.5,90,3);LL,BB=np.meshgrid(np.radians(lgd),np.radians(bgd));W=np.cos(BB)
xyz=np.stack([np.cos(BB)*np.cos(LL),np.cos(BB)*np.sin(LL),np.sin(BB)]).reshape(3,-1)
St=np.array([((np.abs(np.angle(np.exp(1j*(LL-np.radians(c)))))<np.radians(45))*W).ravel() for c in l]);St/=St.sum(1,keepdims=True)
def vec(lo,la):lo,la=np.radians(lo),np.radians(la);return np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
def sky(rings):
  s=np.zeros(xyz.shape[1])
  for (al,ab,th,w,A) in rings:
    t=np.degrees(np.arccos(np.clip(vec(al,ab)@xyz,-1,1)));s+=A*np.exp(-0.5*((t-th)/w)**2)
  return s
R=json.load(open('A_required_rings_results.json'))
def share(v):P=spec(v);return P[2:M+1]/P[2:M+1].sum()
