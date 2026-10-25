import sys,os,pickle,json,warnings;warnings.filterwarnings('ignore')
G='/workspace/cosmology-course/notes/accel/gr';H=os.path.abspath('..');sys.path[:0]=[G+'/valley',G+'/ltb',G,H]
os.chdir(G+'/valley');import toyext as toy;os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np;from ltbcore import photon;from scipy.optimize import brentq
B=json.load(open('../../point_fit.json'))['best'];f=B['core_frac'];uD=10**B['log10_uD'];l0=B['l0'];bf=B['beta_f']
vp=B['v_pec_kms']/299792.458;lp=B['l_pec'];bp=np.degrees(np.arcsin(B['sin_b_pec']))
u=lambda l,b:np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(l0,30.0);dv=u(264.021,48.253);Pp=dv-(dv@A)*A;Pp/=np.linalg.norm(Pp)
pv=vp*u(lp,bp);V=np.array([bf+pv@A,pv@Pp])   # ring-plane velocity (x = away from CoE); out-of-plane part dropped
AREC=[float(os.environ.get('AR',x)) for x in os.environ.get('ARS','0.001,0.0015,0.002,0.003,0.0045,0.0068,0.01').split(',')];TH=np.arange(0,360,2.0)
def job(fm):
    Dt=uD;Mc=max(fm,1e-12)*Dt**3/2;Mfun=lambda r:Mc+0*r;rm=0.02
    while True:
        try:L=toy.LTBM(Mfun,Dt,Dt,rm,rmax=500.,nr=1600,nt=300);break
        except ValueError: rm*=1.5
    T=L.t1;vo=L.frame(T,Dt,0.0,V);R0,Rp0,Rt0,G0,_,_=L.fields(T,Dt)
    aloc=lambda r,t:(L.Rs.ev(r,t)**2*L.Rs.ev(r,t,dx=1)/r**2)**(1/3)
    out={ar:[] for ar in AREC};why=[]
    for th in np.radians(TH):
        b=th+np.pi;so=photon(L,Dt,0.0,b,T,L.tg[0]*1.05);ts=np.geomspace(T,so.t[-1],4000);y=so.sol(ts);al=aloc(y[0],ts)
        why.append((so.t[-1],float(y[0][-1]),float(al[-1])))
        ko=L.kdu(T,Dt,[np.cos(b)/np.sqrt(G0),np.sin(b)/R0],0.0,vo)
        for ar in AREC:
            k=np.where(al<=ar)[0]
            if len(k)==0 or k[0]==0: out[ar].append((np.nan,)*3);continue
            k=k[0];te=brentq(lambda t:aloc(so.sol(t)[0],t)-ar,ts[k],ts[k-1]);ye=so.sol(te)
            out[ar].append((ko/L.kdu(te,ye[0],ye[2:4],ye[4],[0.0,0.0]),te,float(L.Rs.ev(ye[0],te))))
    return fm,dict(res={k:np.array(v) for k,v in out.items()},why=why,t1=T,rmin=rm,tmin=L.tg[0],V=V)
if __name__=='__main__':
    import multiprocessing as mp;R={}
    with mp.Pool(2) as P:
        for fm,r in P.imap_unordered(job,[f,0.0]):
            R[fm]=r;pickle.dump(dict(R=R,params=dict(f=f,uD=uD,l0=l0,bf=bf,vp=vp,lp=lp,bp=bp,V=V.tolist())),open('bf.pkl','wb'))
            for ar,x in r['res'].items(): print(fm,ar,'missing',int(np.isnan(x[:,0]).sum()),flush=True)
