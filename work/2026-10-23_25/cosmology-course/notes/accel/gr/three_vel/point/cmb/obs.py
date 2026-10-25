import numpy as np,healpy as hp,json
T0=2.7255; u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dv=u(264.021,48.253); beta=3362.08e-6/T0; Pp=dv-(dv@A)*A; Pp/=np.linalg.norm(Pp)
th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),Pp)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.array([I[hp.query_disc(ns,v,np.radians(1.0),nest=True)].mean() for v in n])
full=cmb/T0+beta*(n@dv); full-=full[:-1].mean(); np.save('obs_ring.npy',full)
x=full[:-1];print(dict(rms=np.sqrt(np.mean(x**2)),range=x.max()-x.min(),mean_abs=np.mean(abs(x)),med_abs=np.median(abs(x))))
