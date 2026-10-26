# Compensation error metric: max|Resid + Ftarget|, Ftarget=(1-f) * l=2 part of field map (as in p5_compensate.py), plus l>=3 split.
import sys,json,numpy as np,healpy as hp
z=np.load(sys.argv[1]);f=float(sys.argv[2]) if len(sys.argv)>2 else 0.26
NS=hp.npix2nside(len(z['Fmap']));LMX=3*NS-1
a=hp.map2alm(z['Fmap'],lmax=LMX);ll=hp.Alm.getlm(LMX)[0]
Ft=(1-f)*hp.alm2map(a*(ll==2),NS,lmax=LMX)
E=z['Resid']+Ft;ae=hp.map2alm(E,lmax=LMX)
E2=hp.alm2map(ae*(ll==2),NS,lmax=LMX)
out=dict(max_err=float(np.abs(E).max()),rms_err=float(E.std()),max_err_l2=float(np.abs(E2).max()),max_err_lge3=float(np.abs(E-E2-hp.alm2map(ae*(ll<2),NS,lmax=LMX)).max()),max_Ftarget=float(np.abs(Ft).max()))
print(json.dumps(out));json.dump(out,open(sys.argv[1].replace('.npz','_err.json'),'w'))
