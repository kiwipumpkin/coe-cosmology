# Sedov-Taylor self-similar blast wave, j=3, omega=0 (uniform ambient). Integrates the self-similar Euler ODEs
# (u=a r/t V, rho=rho0 G, c^2=(a r/t)^2 Z, a=2/5, xi=r/R_s) inward from strong-shock jump conditions.
import numpy as np,sys,json
from scipy.integrate import solve_ivp
XMIN=3e-3
def sedov(g,n=3000):
    a=0.4
    def rhs(lam,y):
        V,G,Z=y[0],np.exp(y[1]),np.exp(y[2])
        # unknowns s=[gl=dlnG/dlam, dV, dZ]
        M=np.array([[V-1,1,0],[Z/g,V-1,1/g],[(V-1)*(1-g),0,(V-1)/Z]])
        b=np.array([-3*V,V/a-V*V-2*Z/g,2/a-2*V])
        gl,dV,dZ=np.linalg.solve(M,b);return [dV,gl,dZ/Z]
    y0=[2/(g+1),np.log((g+1)/(g-1)),np.log(2*g*(g-1)/(g+1)**2)]
    lam=np.linspace(0,np.log(XMIN),n)
    s=solve_ivp(rhs,(0,lam[-1]),y0,t_eval=lam,rtol=1e-10,atol=1e-12,method='Radau')
    xi=np.exp(s.t);V,G,Z=s.y[0],np.exp(s.y[1]),np.exp(s.y[2])
    I=np.trapezoid((xi**4*(0.5*G*V**2+G*Z/(g*(g-1))))[::-1],xi[::-1])
    xi0=(4*np.pi*a*a*I)**-0.2
    return xi,G,V,Z,xi0
if __name__=='__main__':
    out={}
    for g in (5/3,4/3,1.4):
        xi,G,V,Z,xi0=sedov(g);k=np.searchsorted(xi[::-1],0.1)
        out[str(round(g,4))]=dict(xi0=xi0,G_shock=G[0],compression=(g+1)/(g-1),G_at=[float(np.interp(q,xi[::-1],G[::-1])) for q in (0.5,0.7,0.85,0.95,0.99)],
          center_logslope=float(np.polyfit(np.log(xi[(xi<0.03)&(xi>4e-3)]),np.log(G[(xi<0.03)&(xi>4e-3)]),1)[0]),xi_min=float(xi[-1]))
        np.savez('sedov_g%.4f.npz'%g,xi=xi[::-1],G=G[::-1],V=V[::-1],Z=Z[::-1],xi0=xi0)
    print(json.dumps(out,indent=1));json.dump(out,open('sedov_check.json','w'),indent=1)
