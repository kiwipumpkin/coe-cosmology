import json,numpy as np,zc
from scipy.optimize import minimize
from multiprocessing import Pool
RR=0.5  # D>R with R ~ Migkas depth: enforce D>=1 Gpc (well above a few hundred Mpc)
def f(p):
    a,b,xm,lf,D=p
    if a<=0.05 or b<=0.2 or xm<=0.2 or D<1.0 or lf<0.3:return 1e3
    xp=xm/(a/b)**(1/b);o=zc.run(a,b,xp,10**lf,D);return o.get('chi2',1e3)
def go(x0):
    r=minimize(f,x0,method='Nelder-Mead',options=dict(maxiter=120,xatol=1e-3,fatol=1e-3,initial_simplex=None));return dict(x=r.x.tolist(),chi2=float(r.fun))
if __name__=='__main__':
    S=[[0.5,1,0.9,np.log10(30),42],[0.5,4,0.8,np.log10(30),20],[0.5,1,0.6,np.log10(30),20]]
    with Pool(3) as P:R=P.map(go,S)
    json.dump(R,open('refA.json','w'));print(R)
