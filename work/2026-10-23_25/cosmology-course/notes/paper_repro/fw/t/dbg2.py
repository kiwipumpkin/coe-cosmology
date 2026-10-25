import cmb_gr as g,numpy as np,time
import scipy.integrate as si
orig=si.solve_ivp; cnt=[0]; t0=time.time()
def wrap(f,a,*args,**k):
    def ff(s,y):
        cnt[0]+=1; d=f(s,y)
        if cnt[0]<4 or cnt[0]%3000==0: print(cnt[0],round(time.time()-t0,1),s,np.round(y,6),np.round(d,6),flush=True)
        return d
    k['rtol']=1e-7; k['atol']=1e-9
    r=orig(ff,a,*args,**k); print('nfev',r.nfev,r.status,r.message); return r
g.solve_ivp=wrap
for m in (0,1):
    cnt[0]=0; pa=g.geo((90.,m)); e=g.at(pa,1e-3); print('done',m,len(pa['t']),pa['R'][-1],pa['null_end'],e)
