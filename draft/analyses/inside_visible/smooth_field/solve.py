import numpy as np
from sf_core import *
class Solver:
    def __init__(self,g,s):
        self.g=g; self.s=s
        SA=np.array([g.smooth(a,s) for a in g.A]); Sm=g.smooth(g.m/g.m.mean(),s)
        self.SA=SA; self.D=np.r_[SA,Sm[None]*1.0]          # rows: data + mass row (target 0)
        self.K=self.D@self.D.T
    def field(self,lam,kappa=1e3):
        W=np.r_[np.ones(len(yR)),[kappa]]; Dw=self.D*W[:,None]; K=Dw@Dw.T; t=np.r_[yR,0.0]*W
        u=np.linalg.solve(K+lam*np.eye(len(K)),t); w=Dw.T@u
        return self.g.smooth(w,self.s)                         # delta = S w
    def score(self,d):
        y=self.g.A@d; left=np.sqrt(((yR-y)**2).sum()+TAIL**2)
        return dict(left=float(left),rms=float(np.sqrt(np.mean(d**2))),peak=float(np.abs(d).max()),dmin=float(d.min()),dmax=float(d.max()),
                    net_mass_frac=float((self.g.m*d).sum()/self.g.m.sum()))
    def to_target(self,target=1.0):
        lo,hi=-30,15                                          # log10 lambda bisection: larger lambda -> larger leftover
        for _ in range(60):
            mid=(lo+hi)/2; sc=self.score(self.field(10**mid))
            if sc["left"]>target: hi=mid
            else: lo=mid
        d=self.field(10**lo); sc=self.score(d); sc["lam"]=10**lo; return d,sc
