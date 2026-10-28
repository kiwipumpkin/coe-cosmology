# Orange arc: circle in the plane of (observer, CoE direction, AoE), centred on the CoE at distance D, passing through the observer.
# D is set so the arc crosses our R=1 sphere on the AoE axis: cos(angle CoE-to-crossing) = R/(2D).
import numpy as np
def v(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def lb(u):u=u/np.linalg.norm(u);return np.degrees(np.arctan2(u[1],u[0]))%360,np.degrees(np.arcsin(u[2]))
c=v(180,-15);A=v(260,60)
P=A if c@A>0 else -A                      # AoE end on the CoE side (antipode (80,-60) here)
DR=1/(2*(c@P))                            # D/R
e2=P-(P@c)*c;e2/=np.linalg.norm(e2)       # in-plane unit vector perpendicular to c
Cpos=DR*c
ph=np.linspace(-np.pi,np.pi,4001)
arc=Cpos[:,None]+DR*(-np.outer(c,np.cos(ph))+np.outer(e2,np.sin(ph)))   # passes through origin at ph=0
keep=np.linalg.norm(arc,axis=0)<=2.3; arc=arc[:,keep]
P2=2*(P@c)*c-P                             # mirror crossing in the same plane
if __name__=='__main__':print('D/R=%.3f'%DR,'P1',lb(P),'P2',lb(P2),'chord ends',lb(-P),lb(-P2))
