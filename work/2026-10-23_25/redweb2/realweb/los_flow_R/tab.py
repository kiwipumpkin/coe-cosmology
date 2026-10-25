exec(open('/workspace/redweb2/realweb/los_flow/ltb.py').read().split('def rho')[0])
from scipy.integrate import cumulative_trapezoid as ct
trec=3.2312e-8;Ts=np.geomspace(trec,t0,400);xs=np.linspace(0.55,1.5,191)
Sg=np.array([[S(t,x)[0,0] for x in xs] for t in Ts]);St=np.array([[S(t,x,dx=1)[0,0] for x in xs] for t in Ts])
a=np.array([S(t,1.)[0,0] for t in Ts]);eta=ct(1/a,Ts,initial=0);eta=eta[-1]-eta   # conformal lookback, D/c units
np.savez('/workspace/redweb2/realweb/los_flow_R/tab.npz',Ts=Ts,xs=xs,S=Sg,St=St,a=a,eta=eta,t0=t0)
print('eta(trec) D=',eta[0],'Gpc',42*eta[0],'Rdot at x=1 today',St[-1,np.argmin(abs(xs-1))],'at trec',St[0,np.argmin(abs(xs-1))])
