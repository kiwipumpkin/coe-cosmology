# Flow build with observation distance R as an input. Lagrangian labels; source surface at comoving R; lookback via rescaled conformal time.
import os,sys,time,numpy as np,healpy as hp
sys.path.insert(0,'/workspace/redweb2/realweb');import build as B
B.PC=-B.PC   # convention fix: CoE is at Galactic (180,-15); aC=(0,15) is the ANTI-CoE axis (physics.py, run.py)
sys.argv=[sys.argv[0]]
exec(open('/workspace/redweb2/realweb/los/los.py').read().split("if __name__=='__main__':")[0])  # grids,tocoarse,phi,samp,PC_ (now corrected CoE potential)
from scipy.interpolate import RegularGridInterpolator as RGI
O='/workspace/redweb2/realweb/los_flow_R/';T=np.load(O+'tab.npz');L1=np.load('/workspace/redweb2/realweb/los_flow/ltb.npz')
lTs=np.log(T['Ts']);Sf=RGI((lTs,T['xs']),T['S']);Stf=RGI((lTs,T['xs']),T['St'])
ix1=np.argmin(abs(T['xs']-1));a_=T['S'][:,ix1];adot=T['St'][:,ix1];eta=T['eta'];t0=float(T['t0'])
Afun={'grow':lambda t:np.interp(np.log(t),np.log(L1['T']),L1['Agrow']),'fixed':lambda t:np.interp(np.log(t),np.log(L1['T']),L1['Afix'])}
def t_of_chi(ch,R):return np.exp(np.interp(-ch/R*eta[0],-eta,lTs))   # chi/R = eta(t)/eta(t_rec)
RS=[4.0,8.0,14.0];dch=0.025;GPC_S=42*3.0857e25/2.998e8   # D/c in seconds (unused, units handled analytically)
def comp(fq,ff,R):
    ch=np.arange(0,R,dch)+dch/2;tt=t_of_chi(ch,R);te=t_of_chi(np.r_[ch-dch/2,R],R)   # edges
    dt=-np.diff(te)  # D/c, positive (t decreases with chi)
    ai=np.interp(np.log(tt),lTs,a_);Hi=np.interp(np.log(tt),lTs,adot)/ai   # c/D
    out={k:np.zeros(len(V)) for k in ('isw_grow','isw_fixed','mlr_grow','mlr_fixed','mlf_grow','mlf_fixed','rs_cc')};
    for k in ('al_grow','al_fixed','ac'):out[k]=np.zeros((len(V),3))
    vstat=[np.inf,-np.inf];Ag={n:A(te) for n,A in Afun.items()}
    for i0 in range(0,len(ch),40):
        sl=slice(i0,i0+40);c_=ch[sl];n=len(c_);P=(V[:,None,:]*c_[None,:,None]).reshape(-1,3)
        Xc=P-B.PC;r=np.linalg.norm(Xc,axis=1);x=r/42;uh=Xc/r[:,None];tq=np.repeat(tt[sl][None],len(V),0).ravel();lt=np.log(tq)
        ratio=Sf(np.c_[lt,x])/np.interp(lt,lTs,a_);delta=B.PC+42*ratio[:,None]*uh-P        # comoving (observer frame) displacement of bodies labelled P since t
        # bodies now at P were labelled P-delta: evaluate today's field there (first-order inverse map)
        w=samp(Gq,fq,P-delta);inf=np.linalg.norm(P,axis=1)<9.0
        if inf.any():w[inf]+=samp(Gf,ff,(P-delta)[inf])
        w=w.reshape(len(V),n,4);g=w[...,1:]
        # comoving velocity of bodies relative to observer-centred homologous expansion: d(delta)/dt (finite diff in t)
        tq2=tq*1.001;lt2=np.log(tq2);ratio2=Sf(np.c_[lt2,x])/np.interp(lt2,lTs,a_);ddel=(42*(ratio2-ratio)[:,None]*uh)/(tq2-tq)[:,None]   # Gpc per D/c (comoving)
        Rdot=Stf(np.c_[lt,x]);vstat=[min(vstat[0],Rdot.min()),max(vstat[1],Rdot.max())]
        Rdot=Rdot.reshape(len(V),n);uh3=uh.reshape(len(V),n,3);ddel=ddel.reshape(len(V),n,3)
        for nm in ('grow','fixed'):
            A0=Ag[nm][i0:i0+n];A1=Ag[nm][i0+1:i0+n+1];Am=0.5*(A0+A1)
            out['isw_'+nm]+=(2/c**2)*((A0-A1)[None]*w[...,0]).sum(1)          # dA/dt Phi0 dt (A at nearer edge minus farther edge = forward-time change)
            # moving lens, comoving relative motion:  dPhi/dt|_P = -A grad_X Phi0 . d(delta)/dt   (grad in /Gpc comoving, ddel in Gpc/(D/c))
            out['mlr_'+nm]+=(2/c**2)*(-(Am*dt[sl])[None]*(g*ddel).sum(-1)).sum(1)
            # moving lens with full CoE radial velocity (Chris): dPhi/dt = -v.grad_phys Phi, v=Rdot*c u_hat, grad_phys=grad_X/(a*?)
            # in D/c time units: -(A/a)*Rdot*(u.grad_X Phi0)*42  [grad_X per Gpc, v in Gpc/(D/c) = 42*Rdot]
            out['mlf_'+nm]+=(2/c**2)*(-(Am/ai[sl]*dt[sl])[None]*42*Rdot*(g*uh3).sum(-1)).sum(1)
            lev=(1-c_/R)[None,:,None];gp=g-(g*V[:,None,:]).sum(-1,keepdims=True)*V[:,None,:]
            out['al_'+nm]+=(-2/c**2*lev*dch*Am[None,:,None]*gp).sum(1)
        cC=samp(Gq,PC_,P).reshape(len(V),n,4)[...,1:];gpc=cC-(cC*V[:,None,:]).sum(-1,keepdims=True)*V[:,None,:]
        out['ac']+=(-2/c**2*(1-c_/R)[None,:,None]*dch*gpc).sum(1)
        Hk=Hi[sl]*c/42   # km/s/Gpc
        out['rs_cc']+=(((cC**2).sum(-1))*(2/c**3)*(2/(3*ai[sl]*Hk))*dch).sum(1)
    # SW at source with positions at recombination
    Ps=R*V;Xc=Ps-B.PC;r=np.linalg.norm(Xc,axis=1);uh=Xc/r[:,None];lt=np.full(len(V),np.log(t_of_chi(R-1e-9,R)))
    d=B.PC+42*(Sf(np.c_[lt,r/42])/np.interp(lt,lTs,a_))[:,None]*uh-Ps;w=samp(Gq,fq,Ps-d)[:,0]
    for nm in ('grow','fixed'):out['sw_'+nm]=Afun[nm](np.exp(lt[0]))*w/(3*c**2)
    # source Doppler from anisotropic CoE flow (no web; reported only): -n.(v_src - v_obs - H sep)
    ts=np.exp(lt[0]);H=np.interp(lt[0],lTs,adot)/np.interp(lt[0],lTs,a_);Ss=Sf(np.c_[lt,r/42]);Sd=Stf(np.c_[lt,r/42]);u0=-B.PC/42;s1=np.interp(lt[0],lTs,a_);v1=np.interp(lt[0],lTs,adot)
    va=Sd[:,None]*uh-v1*u0-H*(Ss[:,None]*uh-s1*u0);out['dop_src']=-(va*V).sum(1)
    out['vr_min'],out['vr_max']=vstat;out['delta_src_max']=float(np.linalg.norm(d,axis=1).max())
    return out
if __name__=='__main__':
    ids=sorted(int(f[1:-4]) for f in os.listdir('/workspace/redweb2/realweb/los/comp'));os.makedirs(O+'comp',exist_ok=True)
    for k in ids:
        f=O+f'comp/u{k}.npz'
        if os.path.exists(f):continue
        t=time.time();Rf,Rq=grids(500000+k);Rq=Rq+tocoarse(Rf)
        fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
        res={}
        for R in RS:
            o=comp(fq,ff,R)
            for q,v in o.items():res[f'R{R:g}_{q}']=v
        np.savez(f,**res);print(k,'%.1fs'%(time.time()-t),{q:(float(np.std(v)) if np.ndim(v) else float(v)) for q,v in res.items() if 'R14' in q},flush=True)
