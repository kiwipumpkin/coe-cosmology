# Check: Migkas A/A_all model with Galilean (static-space) velocity addition instead of special-relativistic, at the saved best fit. No refit.
import os,sys,json,numpy as np
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/point')
src=open('point.py').read().split('# multistart')[0];exec(src)
import boostgrad
from scipy.stats import chi2 as C2,norm
names=['core_frac','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec']
J=json.load(open('point_fit.json'));xb=np.array([J['best'][n] for n in names])
wm=np.average(Y,weights=1/S**2);cf=float(np.sum(((Y-wm)/S)**2))
def st(x):
    m=model(x);c=float(np.sum(((m-Y)/S)**2));R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2));sg=float(norm.isf(C2.sf(cf-c,6)/2));return c,R2,sg
print('SR (saved model):',st(xb),'speed',10**xb[1],'mass',xb[0])
def rest_grad_gal(Hpar,Hperp,bf,D,eps=1e-6):
    def w(x):
        r=np.array([x[0]+D,x[1]]);rh=r/np.linalg.norm(r);return np.array([Hpar*x[0],Hperp*x[1]])+bf*rh
    G=np.zeros((2,2))
    for j in range(2):
        e=np.zeros(2);e[j]=eps;G[:,j]=(w(e)-w(-e))/(2*eps)-0  # Galilean: velocity relative to our frame = w - bf*xhat, constant offset drops from gradient
    return G
globals()['rest_grad']=rest_grad_gal
print('Galilean:',st(xb))
pr=json.load(open('prof_c4.json'))
for k in ('4','8','30'):
    x=np.array(pr[k]['x']);print(k,'c Galilean',st(x))
print('flat chi2',cf,'n',len(Y))
