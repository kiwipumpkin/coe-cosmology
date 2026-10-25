import sys,os,json,numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq,least_squares
sys.path.insert(0,'..');from boostgrad import rest_grad
B="/workspace/cosmology/papers/spreadsheet_model/bubble";sys.path.insert(0,B);c=os.getcwd();os.chdir(B)
from bubble_core import IMIG,MA,MLO,MHI;os.chdir(c)
L=IMIG.astype(float);Y=MA;S=0.5*(MLO+MHI);OS=0.049
zz=np.linspace(0.02,0.3,2000);WI=np.trapezoid(zz,zz)/np.trapezoid(zz**2,zz);BP=0.00123
def age(R,M,E): return quad(lambda r:1/np.sqrt(2*M/r+2*E),0,R,limit=200)[0]
def ratio(f,uD):  # units D=1, G=c=1 ; common bang time; returns H_par/H_perp at D
    xin=max(0.,1-1/uD)                # visible-sphere radius = c/H_perp = D/uD
    Menc=lambda R: 0.5*uD**2*(f+OS*max(0.,R**3-xin**3))
    E0=0.5*uD**2-Menc(1.); t0=age(1.,Menc(1.),E0)
    def u(R):
        M=Menc(R); Emin=-M/R*0.999999
        E=brentq(lambda E:age(R,M,E)-t0,Emin,50*uD**2+10)
        return np.sqrt(2*M/R+2*E)
    h=1e-3; return (u(1+h)-u(1-h))/(2*h)/uD
def pred(rat,uD,l0,bf,lp,bp):
    G=rest_grad(rat*uD,uD,bf,1.0); Hpa,Hpe=G[0,0],G[1,1]   # shared-frame rate change (exact SR)
    t=np.radians(L-l0); H=Hpe+(Hpa-Hpe)*np.cos(t)**2
    tt=np.radians(np.arange(360)-l0); Hb=np.mean(Hpe+(Hpa-Hpe)*np.cos(tt)**2)
    h=H/Hb*(1-bf*WI*np.cos(t-np.pi))-BP*WI*np.cos(np.radians(bp))*np.cos(np.radians(L-lp)) # frame dipole: CoE at l0+180, frame moves away from CoE
    return h**2
fs=np.linspace(0.0,0.6,31); us=np.array([0.5,0.8,1,1.5,2,3,5,10.5,20])
RT=np.array([[ratio(f,u) for u in us] for f in fs]); np.save('ratio_grid.npy',RT)
CH=np.full(RT.shape,np.nan);P={}
for i,f in enumerate(fs):
  for j,u in enumerate(us):
    best=None
    for lp0 in [264.,84.]:
      r=least_squares(lambda p:(pred(RT[i,j],u,p[0],p[1],p[2],p[3])-Y)/S,[18,0.004,lp0,48],bounds=([-60,-0.9,0,-90],[100,0.9,360,90]))
      if best is None or r.cost<best.cost: best=r
    CH[i,j]=2*best.cost;P[(i,j)]=best.x
np.save('chi_grid.npy',CH)
i,j=np.unravel_index(np.nanargmin(CH),CH.shape);p=P[(i,j)];m=pred(RT[i,j],us[j],*p)
R2=1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)
cmin=CH.min(); prof_f=np.nanmin(CH,axis=1); ok=fs[prof_f<cmin+1]
out=dict(f=fs[i],uD=us[j],ratio=RT[i,j],l0=p[0],beta_f=p[1],lp=p[2],bp=p[3],chi2=cmin,R2=R2,f_1sig=[ok.min(),ok.max()],
 prof_u=dict(zip(map(str,us),np.nanmin(CH,axis=0))),ratio_vs_f_u1=dict(zip(map(str,fs.round(3)),RT[:,2])))
r=least_squares(lambda q:(pred(RT[i,j],us[j],*q)-Y)/S,p);cov=np.linalg.inv(r.jac.T@r.jac);out['err_l0_bf_lp_bp']=list(np.sqrt(np.diag(cov)))
json.dump(out,open('phys_fit.json','w'),indent=1,default=float);print(json.dumps(out,indent=1,default=float))
import pickle;pickle.dump(dict(fs=fs,us=us,RT=RT,CH=CH,P=P),open('grid.pkl','wb'))
