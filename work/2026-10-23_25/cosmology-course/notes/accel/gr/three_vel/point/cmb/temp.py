# Temperature recombination: shell recombines when its local scale a_loc=(R^2 R'/r^2)^(1/3) (=1 today at r=R) falls to a_rec; T_shell ∝ 1/a_loc.
import sys,os,pickle,warnings;warnings.filterwarnings('ignore')
G='/workspace/cosmology-course/notes/accel/gr';sys.path[:0]=[G+'/valley',G+'/ltb',G]
_h=os.path.dirname(os.path.abspath(__file__));os.chdir(G+'/valley');sys.path.insert(0,_h);import toyext as toy;os.chdir(_h)
import numpy as np;from ltbcore import photon;from scipy.optimize import brentq
F=0.0178;AREC=tuple(np.round(np.geomspace(0.003,0.6,14),5));TH=np.arange(0,360,2.0)
def job(a):
    uD,f=a;Dt=uD;Mc=max(f,1e-12)*Dt**3/2;Mfun=lambda r:Mc+0*r;rm=0.02
    while True:
        try:L=toy.LTBM(Mfun,Dt,Dt,rm,rmax=max(200.,Dt+120),nr=1200,nt=300);break
        except ValueError: rm*=1.5
        if rm>0.5*Dt: return a,'no grid'
    T=L.t1;vo=L.frame(T,Dt,0.0,[0.00123,0.0]);R0,Rp0,Rt0,G0,_,_=L.fields(T,Dt)
    aloc=lambda r,t:(L.Rs.ev(r,t)**2*L.Rs.ev(r,t,dx=1)/r**2)**(1/3)
    res={ar:[] for ar in AREC};tes={ar:[] for ar in AREC};Rem={ar:[] for ar in AREC}
    for th in np.radians(TH):
        b=th+np.pi;so=photon(L,Dt,0.0,b,T,L.tg[0]*1.05);ts=np.geomspace(T,so.t[-1],3000);y=so.sol(ts);al=aloc(y[0],ts)
        vph=[np.cos(b)/np.sqrt(G0),np.sin(b)/R0];ko=L.kdu(T,Dt,vph,0.0,vo)
        for ar in AREC:
            k=np.where(al<=ar)[0]
            if len(k)==0 or k[0]==0: res[ar].append(np.nan);tes[ar].append(np.nan);Rem[ar].append(np.nan);continue
            k=k[0];te=brentq(lambda t:aloc(so.sol(t)[0],t)-ar,ts[k],ts[k-1]);ye=so.sol(te)
            opz=ko/L.kdu(te,ye[0],ye[2:4],ye[4],[0.0,0.0])
            res[ar].append(1/(opz*ar));tes[ar].append(te);Rem[ar].append(float(L.Rs.ev(ye[0],te)))   # T_obs/T_rec * (1/a_rec)... T_obs = T_rec/(1+z); stored = 1/((1+z) a_rec)
    return a,dict(x={k:np.array(v) for k,v in res.items()},te={k:np.array(v) for k,v in tes.items()},Rem={k:np.array(v) for k,v in Rem.items()},rmin=rm,t1=T,twoM_over_D=2*Mc/Dt)
if __name__=='__main__':
    import multiprocessing as mp
    jobs=[(u,f) for u in (1.05,3.0,10.5) for f in (F,0.0)];out={}
    with mp.Pool(6) as P:
        for a,r in P.imap_unordered(job,jobs):
            out[a]=r;pickle.dump(out,open('temp.pkl','wb'))
            if isinstance(r,str): print(a,r,flush=True);continue
            for ar,x in r['x'].items():
                d=x/np.nanmean(x)-1;print(a,ar,'nan',int(np.isnan(x).sum()),'rms',np.nanstd(d),'range',np.nanmax(d)-np.nanmin(d),'mean (1+z)a_rec',1/np.nanmean(x),flush=True)
