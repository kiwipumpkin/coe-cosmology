import os, sys, json, glob, numpy as np
from scipy.stats import chi2, norm
def nsig(chi,k):
    L=-chi2.logsf(chi**2,k)
    if L<30: return float(norm.isf(np.exp(-L)))
    z=np.sqrt(2*L)
    for _ in range(5): z=np.sqrt(2*L-np.log(2*np.pi*z*z))
    return float(z)
from ps_core import *
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import T_nw, A as AN, h, ns, sigma_TH
# ---------- collect scans
rows=[]
for f in glob.glob("log_dual_*.txt"):
    for l in open(f):
        if l.startswith("{"): rows.append(json.loads(l))
lp={}
for f in glob.glob("thr_*.json"):
    t=json.load(open(f)); lp[t["cfg"]]=t
curves={}
for r in rows: curves.setdefault(r["cfg"],[]).append((r["B"],r["left_dual"],r["left_primal"]))
for c in curves: curves[c]=sorted(curves[c])
thr={}
for c,cv in curves.items():
    B=np.array([x[0] for x in cv]); L=np.array([x[1] for x in cv])
    above=B[L>1.0]; Blo=float(above.max()) if len(above) else None      # dual lower bound >1 here => threshold above this B
    t=lp.get(c,{}); thr[c]=dict(B_lower=t.get("B_below",Blo),B_star_dual=t.get("B_star"),B_1x=t.get("B_1x_primal",t.get("B_star")),left_at_B_1x=t.get("left_primal_1x",t.get("left_primal_at_Bstar")),curve=cv)
H_of={"in05":0.5,"in10":1.0,"in20":2.0,"both05":0.5,"both10":1.0}
def fmap(B,p):
    df=B/p["Delta"]; f0=p["f0"]; R0=f0/(1-f0)
    up=((f0+df)/(1-f0-df))/R0-1 if f0+df<1 else np.inf; dn=((f0-df)/(1-f0+df))/R0-1 if f0-df>0 else -1
    return df,up,dn
# ---------- natural scatter & LCDM per domain
tab=[]
for c in ["in05","in10","in20","both05","both10"]:
    if c not in thr or thr[c]["B_1x"] is None: continue
    B=thr[c]["B_1x"]; H=H_of[c]; Vc=H**3; Req=H*(3/(4*np.pi))**(1/3)
    for pn in PATCH:
        p=pinfo(pn); df,up,dn=fmap(B,p); Np=Vc/p["Vp"]; sf=np.sqrt(p["f0"]*(1-p["f0"])/Np)
        tab.append(dict(cfg=c,H=H,patch=pn,dplus=p["dp"],dminus=p["dm"],f0=p["f0"],ratio0=p["f0"]/(1-p["f0"]),B=B,B_lower=thr[c]["B_lower"],
            df=df,ratio_up_pct=100*up,ratio_dn_pct=100*dn,Np=Np,sigma_f=sf,n_sigma_domain=df/sf,
            sigma_LCDM_domain=float(sigma_TH(Req)),n_sigma_LCDM=B/float(sigma_TH(Req))))
json.dump(dict(thresholds=thr,table=tab,patches={k:pinfo(k) for k in PATCH}),open("results.json","w"),indent=1,default=float)
# ---------- coherent (whole-pattern) significance: covariance of the CMB pattern produced by random fields
def mahal(C,y):
    w,V=np.linalg.eigh(C); k=w>w.max()*1e-12; z=V[:,k].T@y; return float(np.sqrt((z**2/w[k]).sum())),int(k.sum())
sig={}
for c in ["in05","both10"]:
    gs=CONFIGS[c](); A=np.concatenate([g.A for g in gs],1); Vc=np.concatenate([np.full(len(g.P),g.H**3) for g in gs])
    out={}
    for pn in ["extreme 50 Mpc","typical ±0.3, 50 Mpc","extreme 30 Mpc","extreme 100 Mpc"]:
        p=pinfo(pn); s2=p["dp"]*(-p["dm"])*p["Vp"]/Vc      # white-noise variance of domain contrast from random patch placement
        C=(A*s2)@A.T; chi,k=mahal(C,yR); rms=float(np.sqrt(np.trace(C)))
        out[pn]=dict(cmb_rms_random_over_floor=rms,needed_over_floor=float(np.linalg.norm(yR)),chi=chi,dof=k,
                     n_sigma=nsig(chi,k))
    # LCDM correlated field on the inner grid (same as smooth_field/checks.py), covariance from realizations
    g=gs[0]; n=g.n; H=g.H; k1=2*np.pi*np.fft.fftfreq(n,d=H)
    KX,KY,KZ=np.meshgrid(k1,k1,np.fft.rfftfreq(n,d=H)*2*np.pi,indexing="ij"); K=np.sqrt(KX**2+KY**2+KZ**2)
    kh=K/1000/h; P=np.where(K>0,AN*kh**ns*T_nw(np.maximum(kh,1e-8))**2,0)/h**3*1e-9
    rng=np.random.default_rng(7); Y=[]
    for i in range(300):
        f=np.fft.irfftn(np.fft.rfftn(rng.normal(size=(n,n,n)))*np.sqrt(P/H**3),s=(n,n,n)); Y.append(g.A@f[g.mask])
    Y=np.array(Y); C=np.cov(Y.T); chi,k=mahal(C,yR)
    out["LCDM"]=dict(cmb_rms_random_over_floor=float(np.sqrt(np.trace(C))),needed_over_floor=float(np.linalg.norm(yR)),chi=chi,dof=k,
                     n_sigma=nsig(chi,k),note="inner region only, 300 realizations")
    sig[c]=out; print(c,json.dumps(out,indent=0),flush=True)
json.dump(sig,open("coherent_significance.json","w"),indent=1)
