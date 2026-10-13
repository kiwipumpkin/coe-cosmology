"""JOINT inside+outside fit with continuity across our 14.3 Gpc boundary and cosmic-web nearest-neighbour rules.
Kernel and geometry reused read-only:
  pattern  = ss_core.unit_field (Newtonian ring-minus-us pull, same as every earlier run), net mass removed with the
             inside_core.FLOC 'local' convention (same as all-open / core-open all-open).
  reach    = core_open.co_core.dmax_dir: bodies anywhere from 0.5 Gpc from us out to max(43 Gpc, far side of the CoE core).
             CoE interior fully open: no wall, no density cap.
Bodies: C = point mass 1e49..2e52 kg (core radius at delta=100); V = void radius 0.5..4 Gpc, depth delta = -0.7*q,
        q in [0.02,1] (q is new: lets the density-trend rule make voids shallower; without that rule the fit may still use it).
Hard rules (all levels): no body sphere straddles our 14.3 Gpc sphere (wholly inside or wholly outside, as all-open);
        mass cores never overlap anything; at level 'free'/'cont' voids also never overlap (as before);
        at rule levels void-void overlap is a SOFT rule (<= 50 % strict).
Soft rules (penalty vectors, dimensionless, 0 = satisfied), weight lam (tunable):
  continuity  shells A = 10-14.3 Gpc, B = 14.3-19 Gpc from us: body number density, void share, mean log size per kind,
              net density contrast must match within tolerances.
  network     wall position of masses, void overlap, void contact needs a mass nearby, void neighbour count 5-7,
              mass links >=3, nearest-mass spacing 1-2 mean void diameters, local mass balance (void deficit vs
              neighbouring masses; neighbourhood net ~0).
  trend       void |delta| <= observed max underdensity at its radius; mass excess <= observed max overdensity times
              rho V over its own domain (radius = half the nearest-neighbour distance)  (density_vs_scale/fits_v2.json,
              chart lines over_noDef / under_noDisp), within tol_dex.
Free (no penalty): orientation of the network vs the CoE axis; overall drift from the cosmic mean density."""
import os,sys,json,numpy as np
ME=os.path.dirname(os.path.abspath(__file__)); CS=os.path.dirname(ME); SS=os.path.join(CS,"solution_space")
os.environ.setdefault("AO","open")
for p in (SS,os.path.join(SS,"all_open"),os.path.join(SS,"asym"),os.path.join(CS,"core_open"),os.path.join(CS,"inside_visible")):
    sys.path.insert(0,p)
_c=os.getcwd()
from ss_core import R,FLOOR,KG,RV,rms,C_G,A_IN_G,GPC,RHO,rho_gpc,TH,unit_field,amps,vol,DELTA_C,GAP,MLO,MHI,RLO,RHI,AXIS
from inside_core import FLOC
from co_core import dmax_dir,core_depth
import asym_core as AC
os.chdir(_c)
GAP_US=0.5; GAP_RING=0.1; QLO=0.02
SC=1/(FLOOR*np.sqrt(360))
F2=json.load(open(os.path.join(CS,"density_vs_scale","fits_v2.json")))
OV=F2["over_noDef"]; UN=F2["under_noDisp"]
def over_trend(r_gpc): return 10**OV["logA"]*(np.asarray(r_gpc)*1000.0)**(-OV["beta"])
def under_trend(r_gpc): return np.minimum(0.95,10**UN["logA"]*(np.asarray(r_gpc)*1000.0)**(-UN["beta"]))
CU=C_G/np.linalg.norm(C_G)
SHA=(10.0,RV); SHB=(RV,19.0)
VA=vol(SHA[1])-vol(SHA[0]); VB=vol(SHB[1])-vol(SHB[0])
RW=(2*RV**3)**(1/3)     # 18.02 Gpc: shell 14.3-18.0 has the same volume as our whole sphere
# continuity zone pairs (A inside, B outside): boundary shells 10-14.3 vs 14.3-19 (the e.g. in the brief), and the
# whole visible sphere 0-14.3 vs the equal-volume shell 14.3-18.0 (added: with only the first pair the fit could tuck
# every body at 5-8 Gpc and pass the rule trivially).
# last column: share of the smooth compensating deficit (net body mass, spread over our visible sphere, FLOC
# convention) that falls in zone A; zone B is outside our sphere and gets none.
ZONES=[(SHA[0],SHA[1],SHB[0],SHB[1],VA,VB,VA/vol(RV)),(0.0,RV,RV,RW,vol(RV),vol(RW)-vol(RV),1.0)]

LEVELS={
 # name: dict(cont=..., net=..., trend=..., lam=...)
 "free":  dict(cont=None,net=None,trend=None,lam=0),
 "cont":  dict(cont=dict(n=1.5,vfrac=0.25,size=0.3,dd=0.005),net=None,trend=None,lam=10),
 "loose": dict(cont=dict(n=2.5,vfrac=0.4,size=0.6,dd=0.02),
               net=dict(wall1=0.75,wall2=1.5,ovmax=0.75,contact=2.0,nbgap=0.6,vnb=(5,7),mnb=3,sp=(0.5,3.0),bal=1.0,net=0.6),
               trend=dict(tol_dex=2.0),lam=3),
 "web":   dict(cont=dict(n=1.5,vfrac=0.25,size=0.3,dd=0.005),
               net=dict(wall1=0.25,wall2=0.5,ovmax=0.5,contact=1.0,nbgap=0.3,vnb=(5,7),mnb=3,sp=(1.0,2.0),bal=0.3,net=0.3),
               trend=None,lam=30),
 "strict":dict(cont=dict(n=1.5,vfrac=0.25,size=0.3,dd=0.005),
               net=dict(wall1=0.25,wall2=0.5,ovmax=0.5,contact=1.0,nbgap=0.3,vnb=(5,7),mnb=3,sp=(1.0,2.0),bal=0.3,net=0.3),
               trend=dict(tol_dex=0.29),lam=30),
}
def sig(x): return 0.5*(1+np.tanh(0.5*np.clip(x,-60,60)))
def _dirs(lon,lat): return np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]

class JN:
    def __init__(self,kinds,level="free",mirror=True,lam=None):
        self.k=np.array(list(kinds)); self.n=len(kinds); self.isC=self.k=="C"; self.isV=~self.isC
        self.L=LEVELS[level]; self.level=level; self.mirror=mirror
        self.lam=self.L["lam"] if lam is None else lam
        lo=np.tile([0,-4*np.pi,-np.pi/2,0,QLO],(self.n,1)); hi=np.tile([1,4*np.pi,np.pi/2,1,1.0],(self.n,1))
        lo[self.isC,4]=0.999999; hi[self.isC,4]=1.0
        self.lo=lo.ravel(); self.hi=hi.ravel()
    def decode(self,p):
        p=p.reshape(self.n,5); u,lon,lat,s,q=p.T; isC=self.isC
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        U=_dirs(lon,lat); dmin=Rc+GAP_US; dmx=dmax_dir(U,Rc); X=U*(dmin+u*(dmx-dmin))[:,None]
        dl=np.where(isC,0.0,-0.7*q)
        M[~isC]=dl[~isC]*rho_gpc(X[~isC])*vol(Rv[~isC])*GPC**3
        return X,M,Rv,Rc,dl
    def pattern_X(self,X,M): return M@unit_field(X)-M.sum()*FLOC
    def pattern(self,p):
        X,M,Rv,Rc,dl=self.decode(p); return self.pattern_X(X,M)
    def leftover(self,p): return rms(R-self.pattern(p))/FLOOR
    # ---------------- hard ----------------
    def hard(self,p,dec=None):
        X,M,Rv,Rc,dl=dec or self.decode(p); d=np.linalg.norm(X,axis=1)
        v=[np.minimum(np.maximum(0,d+Rc-(RV-GAP_RING)),np.maximum(0,RV+GAP_RING-(d-Rc)))]
        if self.n>1:
            i,j=np.triu_indices(self.n,1); dd=np.linalg.norm(X[i]-X[j],axis=1)
            both_v=self.isV[i]&self.isV[j]
            mask=~both_v if self.L["net"] is not None else np.ones_like(both_v)
            v.append(np.where(mask,np.maximum(0,Rc[i]+Rc[j]+GAP-dd),0.0))
        return np.concatenate(v)
    # ---------------- soft ----------------
    def cont_pen(self,dec,t):
        return np.concatenate([self._cont1(dec,t,z) for z in ZONES])
    def _cont1(self,dec,t,z):
        za0,za1,zb0,zb1,VA,VB,cf=z
        X,M,Rv,Rc,dl=dec; d=np.linalg.norm(X,axis=1); w=0.3
        mA=(sig((d-za0)/w) if za0>0 else 1.0)*sig((za1-d)/w); mB=sig((d-zb0)/w)*sig((zb1-d)/w)
        nA=mA.sum()+0.5; nB=mB.sum()+0.5
        pen=[max(0,abs(np.log((nA/VA)/(nB/VB)))-np.log(t["n"]))]
        fA=(mA*self.isV).sum()/nA; fB=(mB*self.isV).sum()/nB
        pen.append(max(0,abs(fA-fB)-t["vfrac"])/t["vfrac"])
        for kk,sz in ((self.isC,np.log10(np.maximum(np.abs(M),1))),(self.isV,np.log10(np.maximum(Rv,1e-3))*2)):
            a=(mA*kk).sum(); b=(mB*kk).sum()
            if a>1e-3 and b>1e-3:
                ma=(mA*kk*sz).sum()/a; mb=(mB*kk*sz).sum()/b; wt=min(1.0,a,b)
                pen.append(wt*max(0,abs(ma-mb)-t["size"])/t["size"])
            else: pen.append(0.0)
        dA=((mA*M).sum()-cf*M.sum())/(RHO*VA*GPC**3); dB=(mB*M).sum()/(RHO*VB*GPC**3)
        pen.append(max(0,abs(dA-dB)-t["dd"])/t["dd"])
        return np.array(pen)
    def net_pen(self,dec,t):
        X,M,Rv,Rc,dl=dec; n=self.n; isC=self.isC; isV=self.isV; pen=[]
        D=np.linalg.norm(X[:,None]-X[None],axis=2)+np.eye(n)*1e9
        Dbar=2*Rv[isV].mean() if isV.any() else 4.0
        iC=np.where(isC)[0]; iV=np.where(isV)[0]; nv=len(iV)
        if nv:
            S=(D[np.ix_(iC,iV)]-Rv[iV][None])/Rv[iV][None]          # mass-to-void-surface, in void radii
            Ss=np.sort(S,axis=1)
            pen+=list(np.maximum(0,Ss[:,0]-t["wall1"])/t["wall1"])
            if nv>=2: pen+=list(np.maximum(0,Ss[:,1]-t["wall2"])/t["wall2"])
        if nv>=2:
            a,b=np.triu_indices(nv,1); A=iV[a]; B=iV[b]; dd=D[A,B]; Rmin=np.minimum(Rv[A],Rv[B])
            ov=(Rv[A]+Rv[B]-dd)/Rmin; pen+=list(np.maximum(0,ov-t["ovmax"])/t["ovmax"])
            touch=sig((Rv[A]+Rv[B]+0.1*Rmin-dd)/(0.05*Rmin))
            P=X[A]+(X[B]-X[A])*(Rv[A]/(Rv[A]+Rv[B]))[:,None]
            if len(iC):
                dm=np.linalg.norm(X[iC][None]-P[:,None],axis=2).min(1)
                pen+=list(touch*np.maximum(0,dm/(t["contact"]*Rmin)-1))
            else: pen+=list(touch)
            # void neighbour count
            lim=(Rv[iV][:,None]+Rv[iV][None])*(1+t["nbgap"]); Dv=D[np.ix_(iV,iV)]
            cnt=sig((lim-Dv)/(0.05*lim)).sum(1)
            lo=min(t["vnb"][0],nv-1); hi=t["vnb"][1]
            pen+=list(np.maximum(0,lo-cnt)/max(lo,1)+np.maximum(0,cnt-hi)/hi)
        if len(iC):
            # mass links: voids whose surface within wall2 radii, masses within 2 Dbar
            c=np.zeros(len(iC))
            if nv: c+=sig((t["wall2"]-S)/0.05).sum(1)
            if len(iC)>1: c+=sig((2*Dbar-D[np.ix_(iC,iC)])/(0.05*Dbar)).sum(1)
            lo=min(t["mnb"],n-1); pen+=list(np.maximum(0,lo-c)/lo)
            if len(iC)>1:
                dn=D[np.ix_(iC,iC)].min(1)
                pen+=list(np.maximum(0,np.log(t["sp"][0]*Dbar/dn))+np.maximum(0,np.log(dn/(t["sp"][1]*Dbar))))
        if nv and len(iC):
            # local balance: each mass shared among voids by proximity of void surface
            g=np.exp(-np.maximum(S,0)**2)+1e-12; a=g/g.sum(1,keepdims=True)
            comp=(M[iC][:,None]*a).sum(0); defi=-M[iV]
            pen+=list(np.maximum(0,np.abs(np.log((comp+1e-3*defi)/defi))-np.log(1+t["bal"])))
        Dn=D.copy(); np.fill_diagonal(Dn,0); W=np.exp(-(Dn/(1.5*Dbar))**2)
        f=(W@M)/(W@np.abs(M)); pen+=list(np.maximum(0,np.abs(f)-t["net"])/t["net"])
        return np.array(pen)
    def domain(self,X):
        D=np.linalg.norm(X[:,None]-X[None],axis=2)+np.eye(len(X))*1e9
        return np.clip(0.5*D.min(1),0.2,10.0) if len(X)>1 else np.full(len(X),10.0)
    def trend_excess(self,dec):
        """dex above the observed max-contrast trend, per body"""
        X,M,Rv,Rc,dl=dec; ex=np.zeros(self.n); rd=self.domain(X)
        ex[self.isV]=np.log10(np.abs(dl[self.isV])/under_trend(Rv[self.isV]))
        cap=over_trend(rd[self.isC])*rho_gpc(X[self.isC])*vol(rd[self.isC])*GPC**3
        ex[self.isC]=np.log10(np.abs(M[self.isC])/cap)
        return ex
    def trend_pen(self,dec,t): return np.maximum(0,self.trend_excess(dec)-t["tol_dex"])
    def soft(self,p,dec=None):
        dec=dec or self.decode(p); out=[]
        if self.L["cont"]: out.append(self.cont_pen(dec,self.L["cont"]))
        if self.L["net"]: out.append(self.net_pen(dec,self.L["net"]))
        if self.L["trend"]: out.append(self.trend_pen(dec,self.L["trend"]))
        return np.concatenate(out) if out else np.zeros(0)
    def asym(self,dec):
        if self.mirror: return np.zeros(0)
        X,M,Rv,Rc,dl=dec; return AC.asym_viol(self.k,X,M,Rv)
    def resid(self,p,lam,lamh):
        dec=self.decode(p); w=self.pattern_X(dec[0],dec[1])
        return np.concatenate([(R-w)*SC,lam*self.soft(p,dec),lamh*self.hard(p,dec),lamh*self.asym(dec)])
    def encode(self,X,M,Rv,q=None):
        k=self.k; X=np.array(X,float); M=np.abs(np.array(M,float)); Rv=np.array(Rv,float); p=np.zeros((self.n,5)); isC=self.isC
        d=np.linalg.norm(X,axis=1); p[:,1]=np.arctan2(X[:,1],X[:,0]); p[:,2]=np.arcsin(np.clip(X[:,2]/d,-1,1))
        p[isC,3]=(np.log10(np.clip(M[isC],MLO,MHI))-np.log10(MLO))/np.log10(MHI/MLO); p[~isC,3]=(np.clip(Rv[~isC],RLO,RHI)-RLO)/(RHI-RLO)
        Rc=np.where(isC,(3*np.clip(M,MLO,MHI)/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3),np.clip(Rv,RLO,RHI))
        U=X/d[:,None]; dmin=Rc+GAP_US; dmx=dmax_dir(U,Rc)
        p[:,0]=np.clip((d-dmin)/(dmx-dmin),0,1); p[:,4]=1.0 if q is None else np.where(isC,1.0,q)
        return np.clip(p.ravel(),self.lo,self.hi)

def fit_one(m,p0,lam=None,iters=4,nfev=40):
    from scipy.optimize import least_squares
    lam=m.lam if lam is None else lam; lamh=30.0; p=np.clip(p0,m.lo,m.hi)
    # continuation on the soft-rule weight: lam/10 -> lam (lets the CMB term shape the layout first)
    for it,lf in enumerate([0.1,1.0,1.0,1.0][:iters]):
        f=lambda q: m.resid(q,lam*lf,lamh)
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=nfev*len(p),xtol=1e-8,ftol=1e-8)
        p=r.x
        if it>=1 and m.hard(p).max(initial=0)<=1e-3 and m.asym(m.decode(p)).max(initial=0)<=1e-3: break
        if it>=1: lamh*=10
    return p

def network_props(k,X,M,Rv):
    """link lengths (in mean void diameters), neighbour counts, orientation, inside/outside split"""
    k=np.array(list(k)); X=np.array(X); M=np.array(M); Rv=np.array(Rv); isC=k=="C"; isV=~isC; n=len(k)
    D=np.linalg.norm(X[:,None]-X[None],axis=2)+np.eye(n)*1e9; Dbar=2*Rv[isV].mean() if isV.any() else 4.0
    iC=np.where(isC)[0]; iV=np.where(isV)[0]
    links=[]
    for a in range(n):
        for b in range(a+1,n):
            if isV[a] and isV[b]: ok=D[a,b]<(Rv[a]+Rv[b])*1.3
            elif isC[a] and isC[b]: ok=D[a,b]<2*Dbar
            else:
                v,c=(a,b) if isV[a] else (b,a); ok=D[a,b]-Rv[v]<0.5*Rv[v]
            if ok: links.append((a,b))
    L=np.array([D[a,b] for a,b in links])/Dbar if links else np.zeros(0)
    deg=np.zeros(n,int)
    for a,b in links: deg[a]+=1; deg[b]+=1
    d=np.linalg.norm(X,axis=1); U=X/d[:,None]
    # orientation: mass-weighted quadrupole of |M| positions about us, major axis vs CoE direction; signed dipole
    wq=np.abs(M); Q=(wq[:,None,None]*(3*U[:,:,None]*U[:,None,:]-np.eye(3)[None])).sum(0)/wq.sum()
    ev,evec=np.linalg.eigh(Q); ax=evec[:,np.argmax(np.abs(ev))]
    dip=(M[:,None]*U).sum(0)/np.abs(M).sum()
    inside=d<RV
    return dict(n_links=len(links),link_len_med=float(np.median(L)) if len(L) else None,link_len=L.tolist(),
        deg=deg.tolist(),deg_mean_C=float(deg[isC].mean()) if isC.any() else None,deg_mean_V=float(deg[isV].mean()) if isV.any() else None,
        void_void_deg=[int(sum(1 for a,b in links if (a==i or b==i) and isV[a] and isV[b])) for i in iV],
        quad_axis_cos_CoE=float(abs(ax@CU)),dipole_cos_CoE=float(dip@CU/max(np.linalg.norm(dip),1e-30)),dipole_amp=float(np.linalg.norm(dip)),
        frac_inside=float(inside.mean()),n_inside=int(inside.sum()),n_core=int((core_depth(X)>0).sum()),
        mass_inside=float(M[inside].sum()),mass_outside=float(M[~inside].sum()),
        frac_CoE_hemi=float((U@CU>0).mean()),Dbar=float(Dbar),dist_med=float(np.median(d)),
        delta_visible=float(M[inside].sum()/(RHO*vol(RV)*GPC**3)),links=links)

# ---------------- LCDM rarity of the fitted contrasts ----------------
_LS=np.loadtxt(os.path.join(CS,"density_vs_scale","lcdm_sigma_tophat.txt"))
def sigma_lcdm(r_gpc):
    """linear top-hat sigma at radius r (Gpc); file is R in Mpc (h=0.70 already applied); log-log interpolation,
    extrapolated beyond the table with the last slope"""
    lr=np.log(_LS[:,0]); ls=np.log(_LS[:,1]); x=np.log(np.asarray(r_gpc)*1000.0)
    sl=(ls[-1]-ls[-2])/(lr[-1]-lr[-2])
    return np.exp(np.where(x<=lr[-1],np.interp(x,lr,ls),ls[-1]+sl*(x-lr[-1])))
def lcdm_nsigma(k,X,M,Rv,dl):
    """per body: contrast over its own region (voids: depth over the void; masses: excess over a sphere of half the
    nearest-neighbour distance) divided by the LCDM linear sigma at that radius"""
    from scipy.stats import norm
    k=np.array(list(k)); X=np.array(X); M=np.array(M); Rv=np.array(Rv); dl=np.array(dl); isC=k=="C"; n=len(k)
    D=np.linalg.norm(X[:,None]-X[None],axis=2)+np.eye(n)*1e9; rd=np.clip(0.5*D.min(1),0.2,10.0)
    r=np.where(isC,rd,Rv); delta=np.where(isC,M/(rho_gpc(X)*vol(rd)*GPC**3),dl)
    ns=np.abs(delta)/sigma_lcdm(r); lp=np.log10(np.maximum(2*norm.sf(ns),1e-300))
    return dict(r=r.tolist(),delta=delta.tolist(),nsigma=ns.tolist(),log10p_each=lp.tolist(),log10p_total=float(lp.sum()))
