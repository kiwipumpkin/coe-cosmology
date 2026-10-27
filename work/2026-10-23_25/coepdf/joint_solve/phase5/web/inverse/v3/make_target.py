# Planck target sky for inverse3.py: SMICA (native Nside 2048) -> alm up to l=191, keep l=2..LT (sharp, LT=130 default; no monopole/dipole),
# synthesised at Nside 64 pixel CENTRES (same sampling as the model kernels, so no ud_grade pixel-window suppression). dT/T units (/2.7255 K).
# SMICA 5' beam ignored (b_l >= 0.997 for l<=130). Also checks Dl(target) vs Planck 2018 TT (Pf) at l=2..LT. Fails on NaN.
import sys,json,numpy as np,healpy as hp
LT=int(sys.argv[1]) if len(sys.argv)>1 else 130;NS=64;LMX=3*NS-1;T0=2.7255e6;W='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/'
m=hp.read_map('/workspace/planck/smica.fits')/2.7255
if not np.all(np.isfinite(m)):raise SystemExit('NONFINITE smica')
a=hp.map2alm(m,lmax=LMX,iter=3);del m
f=np.zeros(LMX+1);f[2:LT+1]=1;T=hp.alm2map(hp.almxfl(a,f),NS,lmax=LMX)
if not np.all(np.isfinite(T)):raise SystemExit('NONFINITE T')
np.save(W+f'target_Tobs_l{LT}.npy',T)
l=np.arange(LMX+1);D=l*(l+1)/(2*np.pi)*hp.anafast(T,lmax=LMX)*T0**2;Da=l*(l+1)/(2*np.pi)*hp.alm2cl(a)*T0**2
Pf=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt');P={int(r[0]):r[1] for r in Pf}
np.save(W+f'target_Dl_l{LT}.npy',D);np.save(W+'smica_Dl_alm.npy',Da)
print('rms target',T.std(),'rms nodip',hp.remove_dipole(T).std());print('l, Dl(target ns64 anafast), Dl(smica alm), Planck TT',[(k,round(D[k],1),round(Da[k],1),round(P[k],1)) for k in (2,10,20,50,100,120,128,130)])
