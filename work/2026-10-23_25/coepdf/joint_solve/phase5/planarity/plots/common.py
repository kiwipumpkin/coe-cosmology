# shared helpers: exec planarity.py up to the result section (no sims, no json write)
import numpy as np,healpy as hp
exec(open('/workspace/coepdf/joint_solve/phase5/planarity/planarity.py').read().split('R={k:')[0])
COE=(180.,-15.);nco=np.array(hp.ang2vec(*COE,lonlat=True))
def eig(alm,l):
    v=full(alm,l);L=LL[l];M=np.array([[np.real(v.conj()@L[i]@L[j]@v) for j in range(3)] for i in range(3)]);M=(M+M.T)/2
    w,e=np.linalg.eigh(M);return w/np.real(v.conj()@v),e
def lb(n):n=n if n[2]>=0 else -n;lo,la=hp.vec2ang(n,lonlat=True);return float(lo[0]),float(la[0])
def tl_about(alm,ax,l):
    v=full(rot_to(alm,*ax),l);p=np.abs(v)**2;return float((p[0]+p[-1])/p.sum()),float(p[l]/p.sum())
ALM={k:hp.map2alm(m,lmax=LM,iter=3) for k,m in sk.items()}
ORDER=['Planck_target','CoE_bg_only','v3_s1_start','v3_s1_model']
LAB={'Planck_target':'Planck SMICA','CoE_bg_only':'CoE bg only','v3_s1_start':'random start','v3_s1_model':'fitted model (v3 s1)'}
