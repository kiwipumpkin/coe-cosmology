# Static (non-expanding) Euclidean space. Matter: radial ballistic from CoE at t=0, Newtonian point core mass (r'' = -M/r^2), no speed cap.
# Light: straight lines at exactly c in the fixed frame. Units D=1 (our distance from CoE today), c=1, G=1. M = f * 0.5 * u^2 (u = our speed today).
# Emission: each emitter emits when it crosses R_rec = 0.5 (fixed sphere, plasma inside is opaque). Observation today at t0.
# Redshift: classical two-body Doppler f_obs/f_emit = (1 - v_obs.k)/(1 - v_src.k), k = unit direction source -> observer.
import sys,numpy as np
sys.path.insert(0,'..');import base as V
from scipy.integrate import quad
from scipy.optimize import brentq
R=0.5
class Model:
    def __init__(s,u,f):
        s.u,s.f=u,f;s.M=f*0.5*u*u;s.E0=0.5*u*u-s.M;M=s.M
        s.t0=s.tr(1.0,s.E0);s.tc=s.tr(R,s.E0)
        s.tmax=(np.pi/2*R**1.5/np.sqrt(2*M) if M>0 else np.inf)   # latest outward crossing of R (shell whose top is exactly R)
    def tr(s,r,E):  # outward travel time from CoE to radius r for energy E (E = v^2/2 - M/r)
        if s.M==0: return r/np.sqrt(2*E)
        return quad(lambda x:1/np.sqrt(max(2*s.M/x+2*E,1e-300)),0,r,limit=200)[0]
    def emitter(s,te):  # energy, speed at R, of the shell crossing R outward at time te
        if s.M==0: v=R/te;return 0.5*v*v,v
        lo=-s.M/R*(1-1e-12);hi=0.5*(R/te)**2*4+10
        E=brentq(lambda E:s.tr(R,E)-te,lo,hi,xtol=1e-14);return E,np.sqrt(max(2*s.M/R+2*E,0))
    def rpos(s,t,E):
        w=np.sign(E)*np.sqrt(2*abs(E));return V.Rfun(t,w,s.M)
def sky(m,psis):
    """psi = look direction measured from the anti-CoE direction (180 = looking at the CoE)."""
    rows=[]
    for ps in psis:
        p=np.radians(ps);n=np.array([np.cos(p),np.sin(p)]);X0=np.array([1.0,0.0])
        b=n@X0;disc=b*b-(1-R*R)
        if disc<0: rows.append(dict(psi=ps,status='misses sphere'));continue
        sd=-b-np.sqrt(disc)
        if sd<=0: rows.append(dict(psi=ps,status='misses sphere'));continue
        te=m.t0-sd;Pe=X0+sd*n
        if te<=0: rows.append(dict(psi=ps,status='would need emission before t=0',s=sd,te=te));continue
        if te>=m.tmax: rows.append(dict(psi=ps,status='no outgoing emitter at that time',s=sd,te=te));continue
        E,ve=m.emitter(te);k=-n;Ph=Pe/R
        num=1-m.u*(k[0]);den=1-ve*(Ph@k)
        st='ok' if (num>0 and den>0) else ('both negative' if (num<0 and den<0) else ('observer outruns wave (num<=0)' if num<=0 else 'source outruns its own wave toward us (den<=0)'))
        rows.append(dict(psi=ps,status=st,s=sd,te=te,Pe=Pe,ve=ve,E=E,num=num,den=den,opz=den/num,rt=m.rpos(m.t0,E),before_our_crossing=te<m.tc))
    return rows
