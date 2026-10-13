"""Fast (numba) version of the jn_core penalties + an exact-structure Jacobian:
 CMB rows: per-body finite differences batched in one kernel call (bodies enter the pattern independently);
 penalty rows: per-column finite differences evaluated in compiled code. Checked against jn_core (numpy) in selftest()."""
import numpy as np, numba as nb
from jn_core import *
import jn_core as J
from web_core import D_G,L_G
ANG=AC.ANG_TOL; DT=np.log(1+AC.D_TOL); ST=np.log(1+AC.S_TOL); PL=AC.PLANE_TOL; NRM=AC.NRM.copy(); CA=AC.CA; SA=AC.SA
CGX,CGY,CGZ=C_G
OV_A=10**OV["logA"]; OV_B=OV["beta"]; UN_A=10**UN["logA"]; UN_B=UN["beta"]
GPC3=GPC**3

@nb.njit(cache=True)
def _sig(x):
    if x>60: x=60.0
    if x<-60: x=-60.0
    return 0.5*(1+np.tanh(0.5*x))

@nb.njit(cache=True)
def pen(X,M,Rv,Rc,dl,isC,fl,t,out):
    """fl: [cont,net,trend,asym,hardvv]; t: tolerances (see TOLV). writes into out, returns count"""
    n=X.shape[0]; k=0
    d=np.empty(n)
    for i in range(n): d[i]=np.sqrt(X[i,0]**2+X[i,1]**2+X[i,2]**2)
    D=np.empty((n,n))
    for i in range(n):
        for j in range(n):
            D[i,j]=np.sqrt((X[i,0]-X[j,0])**2+(X[i,1]-X[j,1])**2+(X[i,2]-X[j,2])**2) if i!=j else 1e9
    nv=0; nc=0; sRv=0.0
    for i in range(n):
        if isC[i]: nc+=1
        else: nv+=1; sRv+=Rv[i]
    Dbar=2*sRv/nv if nv>0 else 4.0
    # ---- hard (scaled separately by caller: first block) ----
    for i in range(n):
        a=max(0.0,d[i]+Rc[i]-(RVv-0.1)); b=max(0.0,RVv+0.1-(d[i]-Rc[i])); out[k]=min(a,b); k+=1
    for i in range(n):
        for j in range(i+1,n):
            if (not isC[i]) and (not isC[j]) and fl[4]==0: out[k]=0.0
            else: out[k]=max(0.0,Rc[i]+Rc[j]+0.02-D[i,j])
            k+=1
    kh=k
    # ---- asym ----
    if fl[3]:
        Ua=np.empty((n,3)); Mi=np.empty((n,3)); sz=np.empty(n)
        for i in range(n):
            Ua[i,0]=X[i,0]/d[i]; Ua[i,1]=X[i,1]/d[i]; Ua[i,2]=abs(X[i,2])/d[i]
            Mi[i,0]=CA*Ua[i,0]+SA*Ua[i,1]; Mi[i,1]=SA*Ua[i,0]-CA*Ua[i,1]; Mi[i,2]=Ua[i,2]
            sz[i]=abs(M[i]) if isC[i] else Rv[i]
        for i in range(n):
            for j in range(i+1,n):
                if isC[i]!=isC[j]: continue
                ch=np.sqrt((Mi[i,0]-Ua[j,0])**2+(Mi[i,1]-Ua[j,1])**2+(Mi[i,2]-Ua[j,2])**2)
                ang=np.degrees(2*np.arcsin(min(1.0,ch/2)))
                dd=abs(np.log(d[j]/d[i])); sl=abs(np.log(sz[j]/sz[i]))
                s=((ang/ANG)**8+(dd/DT)**8+(sl/ST)**8)**0.125
                out[k]=max(0.0,1.17-s); k+=1
        for i in range(n):
            u=abs(X[i,0]*NRM[0]+X[i,1]*NRM[1]+X[i,2]*NRM[2])/d[i]
            pl=np.degrees(np.arcsin(min(1.0,u)))
            out[k]=max(0.0,(PL*1.05-pl)/PL); k+=1
    # ---- continuity ----
    if fl[0]:
      w=0.3
      for zp in range(2):
        za0,za1,zb0,zb1,VAz,VBz,cf=ZONES[zp]
        Mtot=0.0
        for i in range(n): Mtot+=M[i]
        nA=0.5; nB=0.5; vA=0.0; vB=0.0; dA=0.0; dB=0.0
        a_c=0.0;b_c=0.0;sa_c=0.0;sb_c=0.0;a_v=0.0;b_v=0.0;sa_v=0.0;sb_v=0.0
        for i in range(n):
            mA=_sig((d[i]-za0)/w)*_sig((za1-d[i])/w) if za0>0 else _sig((za1-d[i])/w)
            mB=_sig((d[i]-zb0)/w)*_sig((zb1-d[i])/w)
            nA+=mA; nB+=mB; dA+=mA*M[i]; dB+=mB*M[i]
            if isC[i]:
                sz=np.log10(max(abs(M[i]),1.0)); a_c+=mA; b_c+=mB; sa_c+=mA*sz; sb_c+=mB*sz
            else:
                sz=2*np.log10(max(Rv[i],1e-3)); vA+=mA; vB+=mB; a_v+=mA; b_v+=mB; sa_v+=mA*sz; sb_v+=mB*sz
        out[k]=max(0.0,abs(np.log((nA/VAz)/(nB/VBz)))-np.log(t[0])); k+=1
        out[k]=max(0.0,abs(vA/nA-vB/nB)-t[1])/t[1]; k+=1
        if a_c>1e-3 and b_c>1e-3: out[k]=min(1.0,a_c,b_c)*max(0.0,abs(sa_c/a_c-sb_c/b_c)-t[2])/t[2]
        else: out[k]=0.0
        k+=1
        if a_v>1e-3 and b_v>1e-3: out[k]=min(1.0,a_v,b_v)*max(0.0,abs(sa_v/a_v-sb_v/b_v)-t[2])/t[2]
        else: out[k]=0.0
        k+=1
        out[k]=max(0.0,abs((dA-cf*Mtot)/(RHO_*VAz*GPC3)-dB/(RHO_*VBz*GPC3))-t[3])/t[3]; k+=1
    # ---- network ----
    if fl[1]:
        w1=t[4]; w2=t[5]; ovm=t[6]; ctc=t[7]; nbg=t[8]; vlo=t[9]; vhi=t[10]; mnb=t[11]; sp0=t[12]; sp1=t[13]; bal=t[14]; net=t[15]
        iC=np.empty(nc,np.int64); iV=np.empty(nv,np.int64); a=0; b=0
        for i in range(n):
            if isC[i]: iC[a]=i; a+=1
            else: iV[b]=i; b+=1
        S=np.empty((nc,max(nv,1)))
        for a in range(nc):
            for b in range(nv): S[a,b]=(D[iC[a],iV[b]]-Rv[iV[b]])/Rv[iV[b]]
        if nv>=1:
            for a in range(nc):
                s1=1e9; s2=1e9
                for b in range(nv):
                    v=S[a,b]
                    if v<s1: s2=s1; s1=v
                    elif v<s2: s2=v
                out[k]=max(0.0,s1-w1)/w1; k+=1
                if nv>=2: out[k+nc-1]=max(0.0,s2-w2)/w2
            if nv>=2: k+=nc
        if nv>=2:
            for a in range(nv):
                for b in range(a+1,nv):
                    A=iV[a]; B=iV[b]; Rm=min(Rv[A],Rv[B]); dd=D[A,B]
                    out[k]=max(0.0,(Rv[A]+Rv[B]-dd)/Rm-ovm)/ovm; k+=1
            for a in range(nv):
                for b in range(a+1,nv):
                    A=iV[a]; B=iV[b]; Rm=min(Rv[A],Rv[B]); dd=D[A,B]
                    touch=_sig((Rv[A]+Rv[B]+0.1*Rm-dd)/(0.05*Rm)); f=Rv[A]/(Rv[A]+Rv[B])
                    px=X[A,0]+(X[B,0]-X[A,0])*f; py=X[A,1]+(X[B,1]-X[A,1])*f; pz=X[A,2]+(X[B,2]-X[A,2])*f
                    if nc>0:
                        dm=1e9
                        for c in range(nc):
                            j=iC[c]; e=np.sqrt((X[j,0]-px)**2+(X[j,1]-py)**2+(X[j,2]-pz)**2)
                            if e<dm: dm=e
                        out[k]=touch*max(0.0,dm/(ctc*Rm)-1)
                    else: out[k]=touch
                    k+=1
            lo=min(vlo,nv-1.0)
            for a in range(nv):
                cnt=0.0; A=iV[a]
                for b in range(nv):
                    if b==a: continue
                    B=iV[b]; lim=(Rv[A]+Rv[B])*(1+nbg); cnt+=_sig((lim-D[A,B])/(0.05*lim))
                out[k]=max(0.0,lo-cnt)/max(lo,1.0)+max(0.0,cnt-vhi)/vhi; k+=1
        if nc>=1:
            lo=min(mnb,n-1.0)
            for a in range(nc):
                c=0.0
                for b in range(nv): c+=_sig((w2-S[a,b])/0.05)
                for b in range(nc):
                    if b!=a: c+=_sig((2*Dbar-D[iC[a],iC[b]])/(0.05*Dbar))
                out[k]=max(0.0,lo-c)/lo; k+=1
            if nc>=2:
                for a in range(nc):
                    dn=1e9
                    for b in range(nc):
                        if b!=a and D[iC[a],iC[b]]<dn: dn=D[iC[a],iC[b]]
                    dn=max(dn,1e-6); out[k]=max(0.0,np.log(sp0*Dbar/dn))+max(0.0,np.log(dn/(sp1*Dbar))); k+=1
        if nv>=1 and nc>=1:
            comp=np.zeros(nv)
            for a in range(nc):
                g=np.empty(nv); gs=0.0
                for b in range(nv):
                    s=max(S[a,b],0.0); g[b]=np.exp(-s*s)+1e-12; gs+=g[b]
                for b in range(nv): comp[b]+=M[iC[a]]*g[b]/gs
            for b in range(nv):
                de=-M[iV[b]]; out[k]=max(0.0,abs(np.log((comp[b]+1e-3*de)/de))-np.log(1+bal)); k+=1
        for i in range(n):
            s1=0.0; s2=0.0
            for j in range(n):
                dd=0.0 if j==i else D[i,j]
                ww=np.exp(-(dd/(1.5*Dbar))**2); s1+=ww*M[j]; s2+=ww*abs(M[j])
            out[k]=max(0.0,abs(s1/s2)-net)/net; k+=1
    # ---- density trend ----
    if fl[2]:
        tol=t[16]
        for i in range(n):
            if isC[i]:
                dn=1e9
                for j in range(n):
                    if D[i,j]<dn: dn=D[i,j]
                rd=min(max(0.5*dn,0.2),10.0) if n>1 else 10.0
                rc=np.sqrt((X[i,0]-CGX)**2+(X[i,1]-CGY)**2+(X[i,2]-CGZ)**2); rho=RHO_*np.exp(-(rc-D_G)/L_G)
                cap=OV_A*(rd*1000.0)**(-OV_B)*rho*4.0/3.0*np.pi*rd**3*GPC3
                ex=np.log10(abs(M[i])/cap)
            else:
                ut=min(0.95,UN_A*(Rv[i]*1000.0)**(-UN_B)); ex=np.log10(abs(dl[i])/ut)
            out[k]=max(0.0,ex-tol); k+=1
    return k,kh

RVv=float(RV); VA_=float(J.VA); VB_=float(J.VB); RHO_=float(RHO)
ZONES=np.array(J.ZONES)

def tolvec(L):
    t=np.zeros(17)
    if L["cont"]: c=L["cont"]; t[0:4]=[c["n"],c["vfrac"],c["size"],c["dd"]]
    if L["net"]: q=L["net"]; t[4:16]=[q["wall1"],q["wall2"],q["ovmax"],q["contact"],q["nbgap"],q["vnb"][0],q["vnb"][1],q["mnb"],q["sp"][0],q["sp"][1],q["bal"],q["net"]]
    if L["trend"]: t[16]=L["trend"]["tol_dex"]
    return t

class JF(JN):
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw)
        self.fl=np.array([self.L["cont"] is not None,self.L["net"] is not None,self.L["trend"] is not None,not self.mirror,self.L["net"] is None],np.int64)
        self.t=tolvec(self.L); self.buf=np.zeros(20000)
        X,M,Rv,Rc,dl=self.decode(np.random.default_rng(0).uniform(self.lo,self.hi)); kk,kh=pen(X,M,Rv,Rc,dl,self.isC,self.fl,self.t,self.buf); self.np_=kk; self.nh=kh
    def pens(self,dec):
        X,M,Rv,Rc,dl=dec; out=np.zeros(self.np_); pen(X,M,Rv,Rc,dl,self.isC,self.fl,self.t,out); return out
    def hard(self,p,dec=None): return self.pens(dec or self.decode(p))[:self.nh]
    def softall(self,dec): return self.pens(dec)[self.nh:]
    def wvec(self,lam,lamh):
        """weights per penalty row: hard and asym get lamh, continuity/network/trend get lam"""
        w=np.full(self.np_,float(lam)); w[:self.nh]=lamh
        if not self.mirror:
            na=sum(1 for i in range(self.n) for j in range(i+1,self.n) if self.isC[i]==self.isC[j])+self.n
            w[self.nh:self.nh+na]=lamh
        return w
    def split(self,dec):
        """(hard, asym, soft) arrays"""
        v=self.pens(dec); na=0
        if not self.mirror: na=sum(1 for i in range(self.n) for j in range(i+1,self.n) if self.isC[i]==self.isC[j])+self.n
        return v[:self.nh],v[self.nh:self.nh+na],v[self.nh+na:]
    def fun(self,p,w):
        dec=self.decode(p); return np.nan_to_num(np.concatenate([(R-self.pattern_X(dec[0],dec[1]))*SC,w*self.pens(dec)]),nan=1e3,posinf=1e3,neginf=-1e3)
    def jac(self,p,w):
        n=self.n; P=p.reshape(n,5); h=np.maximum(1e-7,1e-6*(self.hi-self.lo)).reshape(n,5)
        X,M,Rv,Rc,dl=self.decode(p); base=self.pens((X,M,Rv,Rc,dl)); F0=unit_field(X)
        Jm=np.zeros((360,5*n)); Jp=np.zeros((self.np_,5*n)); out=np.zeros(self.np_)
        for t in range(5):
            Q=P.copy(); hh=h[:,t].copy()
            up=Q[:,t]+hh>self.hi.reshape(n,5)[:,t]; hh[up]=-hh[up]; Q[:,t]+=hh
            X2,M2,Rv2,Rc2,dl2=self.decode(Q.ravel()); F2=unit_field(X2)
            dW=(M2[:,None]*F2-M[:,None]*F0)-((M2-M)/1.0)[:,None]*FLOC[None]
            Jm[:,t::5]=(-dW*SC/hh[:,None]).T
            for b in range(n):
                Xb=X.copy(); Mb=M.copy(); Rvb=Rv.copy(); Rcb=Rc.copy(); dlb=dl.copy()
                Xb[b]=X2[b]; Mb[b]=M2[b]; Rvb[b]=Rv2[b]; Rcb[b]=Rc2[b]; dlb[b]=dl2[b]
                pen(Xb,Mb,Rvb,Rcb,dlb,self.isC,self.fl,self.t,out)
                Jp[:,5*b+t]=w*(out-base)/hh[b]
        return np.nan_to_num(np.vstack([Jm,Jp]),nan=0.0,posinf=0.0,neginf=0.0)

def fit(m,p0,lam=None,stages=(0.1,1,1,1),nfev=150,lamh=30.0):
    from scipy.optimize import least_squares
    lam=m.lam if lam is None else lam; p=np.clip(p0,m.lo,m.hi)
    for it,lf in enumerate(stages):
        w=m.wvec(lam*lf,lamh)
        r=least_squares(lambda q:m.fun(q,w),p,jac=lambda q:m.jac(q,w),bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=nfev,xtol=1e-9,ftol=1e-9)
        p=r.x
        if it>=1:
            h,a,_=m.split(m.decode(p))
            if h.max(initial=0)<=1e-3 and a.max(initial=0)<=1e-3 and lf==stages[-1]: break
            if h.max(initial=0)>1e-3 or a.max(initial=0)>1e-3: lamh*=10
    return p

def selftest():
    rng=np.random.default_rng(3)
    for lv in LEVELS:
        for mir in (True,False):
            m=JF("C"*6+"V"*7,lv,mir); m0=JN("C"*6+"V"*7,lv,mir); p=rng.uniform(m.lo,m.hi); dec=m.decode(p)
            h,a,s=m.split(dec); s0=m0.soft(p,dec); h0=m0.hard(p,dec); a0=m0.asym(dec)
            print(lv,mir,len(s),len(s0),np.abs(np.sort(s)-np.sort(s0)).max(initial=0),np.abs(h-h0).max(),np.abs(np.sort(a)-np.sort(a0)).max(initial=0))
    m=JF("C"*5+"V"*5,"web",False); p=rng.uniform(m.lo,m.hi); w=m.wvec(10,30); Ja=m.jac(p,w)
    from scipy.optimize._numdiff import approx_derivative
    Jn=approx_derivative(lambda q:m.fun(q,w),p,method="2-point",rel_step=1e-7)
    print("jac rel err CMB",np.abs(Ja[:360]-Jn[:360]).max()/np.abs(Jn[:360]).max(),"pen",np.abs(Ja[360:]-Jn[360:]).max()/max(1e-12,np.abs(Jn[360:]).max()))
if __name__=="__main__": selftest()
