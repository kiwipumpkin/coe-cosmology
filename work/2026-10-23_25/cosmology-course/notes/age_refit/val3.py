# FLRW check of the 3000 K criterion: homogeneous -> every direction must give 1+z = 3000/2.725 = 1100.92 exactly
import numpy as np,ltb2 as T,forward as F
P=F.params_from('exp',40,0.6,1.0,13,8); L=T.LTB(P,'fwd'); o=F.integrate(P,True); io=o['io']; xo=np.log(F.R0[io]); tn=o['t']*F.C
from scipy.integrate import solve_ivp
ro=L.rho_rel(tn,np.array([xo]))[0]; LN=3*np.log(3000/2.725)
for th in np.radians([0,90,180]):
    e=lambda t,y: np.log(L.rho_rel(t,np.array([y[0]]))[0]/ro)-LN; e.terminal=True
    s=solve_ivp(T.rhs(L,-L.R(tn,np.array([xo]))[0]*np.sin(th)),(tn,1e-6*tn),[xo,0,0,-np.cos(th)],method='LSODA',rtol=1e-11,atol=1e-13,events=e)
    print('theta',np.degrees(th),'1+z',np.exp(s.y_events[0][0][2]),'target',3000/2.725)
