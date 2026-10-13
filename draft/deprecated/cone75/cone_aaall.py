import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("aaall.py").read().split("def plot")[0])   # A(sg,mf,R,D), ang step .25
i180=721; th=np.radians(ang[:i180]); dth=np.gradient(th)
phi=np.linspace(0,np.pi,721)
def kernel(alpha,taper):
    a=np.radians(alpha); K=np.zeros((i180,i180))
    for i,t0 in enumerate(th):
        cs=np.cos(t0)*np.cos(th)[:,None]+np.sin(t0)*np.sin(th)[:,None]*np.cos(phi)[None,:]
        sep=np.arccos(np.clip(cs,-1,1)); inside=sep<=a
        wt=np.where(inside,np.cos(sep/a*np.pi/2)**2 if taper else 1.0,0.0)
        K[i]=np.trapezoid(wt,phi,axis=1)*np.sin(th)*dth
    K/=K.sum(1,keepdims=True); return K
Ks={"raw":None,"60° top-hat":kernel(60,0),"75° top-hat":kernel(75,0),"60° cos taper":kernel(60,1),"75° cos taper":kernel(75,1)}
def sm(y,K):
    h=y[:i180] if K is None else K@y[:i180]; return h
def full(h): return np.concatenate([h,h[::-1][1:]])
P1=np.cos(th); P2=1.5*np.cos(th)**2-.5; ws=np.sin(th)*dth
def leg(h): m=(h*ws).sum()/ws.sum(); return m,1.5*(h*P1*ws).sum()/ws.sum()*2/2,2.5*(h*P2*ws).sum()/ws.sum()
def metr(h):
    m,a1,a2=leg(h); deg=np.degrees(th)
    return "0°=%.3f 90°=%.3f 180°=%.3f | max %.3f@%.0f min %.3f@%.0f | mean %.3f dip %+.3f quad %+.3f"%(h[0],h[360],h[-1],h.max(),deg[h.argmax()],h.min(),deg[h.argmin()],m,a1,a2)
sty={"raw":("-",2.4),"60° top-hat":("-",1.6),"75° top-hat":("--",1.6),"60° cos taper":(":",1.8),"75° cos taper":("-.",1.4)}
for mf in (0.5,1,2,4):
    f,ax=plt.subplots(figsize=(10,5.5)); ax.axhline(1,color="0.5",lw=.8)
    for sg,col,nm in ((1,"#c0392b","expanding"),(-1,"#2e5fa8","collapsing")):
        y=A(sg,mf)
        for k,K in Ks.items():
            h=sm(y,K); ls,lw=sty[k]; ax.plot(ang,full(h),color=col,ls=ls,lw=lw,alpha=1 if k=="raw" else .85,label="%s, %s"%(nm,k))
            print("M=%g %s %-14s %s"%(mf,nm,k,metr(h)))
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
    ax.set_ylabel("A / A_all  (Chris sheet formula)")
    ax.set_title("A/A_all raw vs Migkas-cone averaged, point mass %g× (strength %.3g), D = 58.2, R = 14.3 Gpc\ntop-hat = solid-angle weighted; cos taper = weight cos²(90°·θ₁/θc) (Migkas 2020 Sect. 4.3 σ/cos)"%(mf,0.247*mf),fontsize=10)
    ax.legend(frameon=False,fontsize=8,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13))
    f.savefig("cone_aaall_M%g.png"%mf,dpi=120,bbox_inches="tight"); plt.close(f)
