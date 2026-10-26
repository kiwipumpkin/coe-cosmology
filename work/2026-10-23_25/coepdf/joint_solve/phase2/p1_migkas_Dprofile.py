# Phase2 A1: profile Migkas chi2 vs fixed D (eta, mc, l0 free; CatWISE delta lock d=RV), reusing catwise-eta-density-migkas.py unchanged.
import sys,json,time,numpy as np
import importlib.util as iu
sp=iu.spec_from_file_location('cw','/workspace/cosmology-course/notes/catwise-eta-density-migkas.py');cw=iu.module_from_spec(sp);sp.loader.exec_module(cw)
from scipy.optimize import minimize
OUT='/workspace/coepdf/joint_solve/phase2/'
import os;Ds=[float(os.environ['DFIX'])]
res={}
for D in Ds:
    best=(np.inf,None);t=time.time()
    for u0 in ([np.log10(4.34),19.9,np.log10(0.05)],[np.log10(1.0),15.0,np.log10(1.0)]):
        r=minimize(lambda u:cw.chi2(D,10**u[2],10**u[0],u[1],cw.RV),u0,method='Nelder-Mead',options=dict(maxiter=160,xatol=1e-4,fatol=1e-4,adaptive=True))
        if r.fun<best[0]:best=(r.fun,r.x)
    mc,l0,eta=10**best[1][0],best[1][1],10**best[1][2]
    full,rr,LD,L=cw.model_bins(D,eta,mc,l0,cw.RV);m=full[cw.IMIG]
    res[str(D)]=dict(D=D,chi2=float(best[0]),r2_pearson=float(np.corrcoef(cw.Y,m)[0,1]**2),eta=eta,mc=mc,l0=float(l0),L=float(L),age_Gyr=float(rr['age']),Ht=float(rr['Ht']),curve=full.tolist())
    print(D,{k:v for k,v in res[str(D)].items() if k!='curve'},'%.0fs'%(time.time()-t),flush=True)
    json.dump(res,open(OUT+f'p1_migkas_D{D:g}.json','w'))
