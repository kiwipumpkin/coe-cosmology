import cmb_gr as g,numpy as np,time,dyn
import scipy.integrate as si
orig=si.solve_ivp
cnt=[0];t0=time.time()
def wrap(f,*a,**k):
    def ff(s,y):
        cnt[0]+=1
        if cnt[0]%2000==0: print(cnt[0],round(time.time()-t0,1),s,y,flush=True)
        return f(s,y)
    return orig(ff,*a,**k)
g.solve_ivp=wrap
pa=g.geo((90.,1)); print('done',len(pa['t']),pa['R'][-1],pa['status'])
