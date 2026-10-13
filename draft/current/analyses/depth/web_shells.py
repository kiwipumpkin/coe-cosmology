"""(B) Random-web effect vs distance and (C) analytic falloff.
B1: regenerate the void_mix k=0, N=200, level b, 50%-void layouts (seed 2, identical to ../shell_ensemble/void_mix raw_k0_N200_b_fv50.npz;
    checked) and split each layout's pattern into distance shells.
B2: same object population (uniform in volume, region s>R_vis, CoE distance >= D-d) extended far beyond 42.9 Gpc, unit masses,
    to get the share of the effect inside R for a web that does not stop at 42.9 Gpc.
'Effect' = power added to the leftover (mean square of the web pattern; E[leftover^2] = base^2 + E[web^2], checked earlier),
which adds object by object because signs/positions are random. Also given per mode group."""
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); VM="/workspace/clump_scan/shell_ensemble/void_mix"
vm={"__file__":os.path.join(VM,"ensemble_void_mix.py"),"__name__":"imp"}
exec(open(vm["__file__"]).read().split('\nif __name__')[0],vm); os.chdir(HERE)
field3,amps,rms,RV,FLOOR,KG,LEVELS,radius_gpc,_draw,C,R=[vm[k] for k in ("field3","amps","rms","RV","FLOOR","KG","LEVELS","radius_gpc","_draw","C","g")][:10]+[vm["g"]["R"]]
A_IN=vm["g"]["g"]["A_IN"]
th=vm["g"]["g"]["g"]["th"]
def modepow(W):   # W (...,360) -> power per mode m=1..12 (A_m^2/2, sums to mean square for m<=12)
    m=np.arange(1,13)[:,None]; c=2*np.mean(W[...,None,:]*np.cos(m*th),-1); s=2*np.mean(W[...,None,:]*np.sin(m*th),-1); return 0.5*(c*c+s*s)
GROUPS={"all (leftover RMS)":None,"m=1-2":[0,1],"m=3-5":[2,3,4]}

# ---------------- B1
def layouts(nreal=1000,N=200,lev="b",fv=0.5,seed=2):
    rng=np.random.default_rng(seed); lo,hi=LEVELS[lev]
    for _ in range(nreal):           # identical rng sequence to void_mix realization() for k=0
        nv=int(np.clip(round(fv*N),1,N-1)); nc=N-nv
        Mc=10**rng.uniform(lo,hi,nc); Mv0=10**rng.uniform(lo,hi,nv); Mv=Mv0*Mc.sum()/Mv0.sum(); M=np.r_[Mc,-Mv]
        P=np.empty((N,3)); P[:nc]=_draw(nc,rng,C,np.ones(nc),1,True); P[nc:]=_draw(nv,rng,C,np.ones(nv),-1,True)
        yield M,P
EDGES=[14.3,20,25,30,35,42.9]
def b1():
    raw=np.load(os.path.join(VM,"raw_k0_N200_b_fv50.npz")); base=rms(R)
    nS=len(EDGES)-1; Pw=np.zeros((1000,nS)); Pm=np.zeros((1000,nS,12)); M2=np.zeros((1000,nS)); dleft=np.zeros((1000,nS)); nobj=np.zeros((1000,nS))
    per_obj=[]; check=[]
    for i,(M,P) in enumerate(layouts()):
        s=np.linalg.norm(P,axis=1); f=(M/KG)[:,None]*field3(np.arctan2(P[:,1],P[:,0]),np.arcsin(P[:,2]/s),s)
        check.append(abs(rms(R-f.sum(0))-raw["rms"][i]))
        sh=np.digitize(s*RV,EDGES)-1
        for j in range(nS):
            w=f[sh==j].sum(0); Pw[i,j]=np.mean(w*w); Pm[i,j]=modepow(w); M2[i,j]=np.sum(M[sh==j]**2); nobj[i,j]=np.sum(sh==j)
            dleft[i,j]=rms(R-w)-base
        po=np.mean(f*f,1); pm=modepow(f)
        per_obj.append(np.c_[s*RV,po,pm])
    per_obj=np.concatenate(per_obj)
    out=dict(max_reproduction_error=float(max(check)),edges_gpc=EDGES,shells=[])
    for j in range(nS):
        unit=np.sqrt(Pw[:,j].sum()/M2[:,j].sum())*1e51/FLOOR      # RMS pattern per 1e51 kg of object mass, x floor
        unit_m=np.sqrt(Pm[:,:,:][:,j].sum(0)/M2[:,j].sum())*1e51/FLOOR
        out["shells"].append(dict(range_gpc=[EDGES[j],EDGES[j+1]],mean_objects=float(nobj[:,j].mean()),
            share_of_power=float(Pw[:,j].sum()/Pw.sum()),share_m12=float(Pm[:,j,0:2].sum()/Pm[:,:,0:2].sum()),share_m35=float(Pm[:,j,2:5].sum()/Pm[:,:,2:5].sum()),
            rms_per_1e51kg_x_floor=float(unit),mode_rms_per_1e51kg_x_floor=(unit_m*np.sqrt(2)).tolist(),   # sqrt(2): amplitude A_m
            median_abs_change_in_leftover_x_floor=float(np.median(np.abs(dleft[:,j]))/FLOOR),mean_change_in_leftover_x_floor=float(dleft[:,j].mean()/FLOOR)))
    out["sum_of_shell_powers_over_total_power"]=float(Pw.sum()/np.sum([np.mean((0)) for _ in [0]]) if False else 0)
    return out,per_obj

# ---------------- B2: unit masses, uniform in volume, region s>s0, CoE distance >= D-d, out to S_BIG
def b2(s0_gpc,S_BIG=400.0,n=600000,seed=7):
    rng=np.random.default_rng(seed); s0=s0_gpc/RV; S=S_BIG/RV
    s=np.exp(rng.uniform(np.log(s0),np.log(S),n)); wvol=s**3            # log-uniform radius, weight = volume per dln s
    u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1)[:,None]; P=s[:,None]*u
    ok=np.linalg.norm(P-C,axis=1)>=A_IN; s,P,wvol=s[ok],P[ok],wvol[ok]
    po=np.empty(len(s)); pm=np.empty((len(s),12))
    for i in range(0,len(s),20000):
        sl=slice(i,i+20000); f=field3(np.arctan2(P[sl,1],P[sl,0]),np.arcsin(P[sl,2]/s[sl]),s[sl]); po[sl]=np.mean(f*f,1); pm[sl]=modepow(f)
    return np.c_[s*RV,wvol*po,wvol[:,None]*pm]
def cum_depth(tab,cols,levels=(0.9,0.95)):
    d=tab[:,0]; w=tab[:,cols].sum(1) if np.ndim(cols) else tab[:,cols]; o=np.argsort(d); c=np.cumsum(w[o])/w.sum()
    return {f"R{int(l*100)}_gpc":float(d[o][np.searchsorted(c,l)]) for l in levels},(d[o],c)
COLS={"all (leftover RMS)":1,"m=1-2":[2,3],"m=3-5":[4,5,6]}

if __name__=="__main__":
    res={}
    o1,per_obj=b1(); o1.pop("sum_of_shell_powers_over_total_power",None); res["B1_ensemble_k0_N200_b_fv50"]=o1
    res["B1_cumulative_depth_within_42.9"]={g:cum_depth(per_obj,c)[0] for g,c in COLS.items()}
    curves={}
    res["B2_extended_web"]={}
    for s0 in (14.3,15.0,15.5):
        tab=b2(s0); res["B2_extended_web"][f"inner_edge_{s0}"]={g:cum_depth(tab,c)[0] for g,c in COLS.items()}
        if s0==15.0:
            for g,c in COLS.items():
                d,cc=cum_depth(tab,c)[1]; idx=np.searchsorted(d,np.linspace(15.0,120,300)); idx=np.clip(idx,0,len(d)-1); curves[g]=dict(R=np.linspace(15.0,120,300).tolist(),share=cc[idx].tolist())
            # per-unit-mass falloff by mode (direction averaged), for (C)
            bins=np.array([15,17,20,25,30,36,43,50,60,80,100,150,200])
            fall=[]
            for a,b in zip(bins[:-1],bins[1:]):
                sel=(tab[:,0]>=a)&(tab[:,0]<b); wv=(tab[sel,0]/RV)**3
                fall.append(dict(r_mid_gpc=float(np.sqrt(a*b)),mode_amp=(np.sqrt(2*tab[sel][:,2:14].sum(0)/wv.sum())).tolist()))
            res["C_falloff_unit_k"]=fall
    # (C) analytic: amplitude per unit mass ~ (R_vis/r)^(m+1) for m>=2, (R_vis/r)^4 for m=1 (l=1 uniform pull removed);
    # uniform population from r0 outward: share inside R = 1-(r0/R)^(2p-3) with p the falloff power -> m=1: 5, m>=2: 2m-1
    an={}
    for m in range(1,6):
        p=4 if m==1 else m+1; e=2*p-3
        an[f"m={m}"]=dict(amp_falloff_power=p,share_inside_R=f"1-(r0/R)^{e}",R90_over_r0=10**(1/e),R95_over_r0=20**(1/e),
                         R90_gpc_r0_14p3=14.3*10**(1/e),R95_gpc_r0_14p3=14.3*20**(1/e))
    res["C_analytic"]=an
    # measured large-r slopes from B2 falloff (log-log between 60-200 Gpc)
    f=res["C_falloff_unit_k"]; r=np.array([x["r_mid_gpc"] for x in f]); A=np.array([x["mode_amp"] for x in f])
    sel=(r>55)&(r<180); res["C_measured_power_law_slopes_60_180gpc"]={f"m={m+1}":float(np.polyfit(np.log(r[sel]),np.log(A[sel,m]),1)[0]) for m in range(6)}
    sel=(r>15)&(r<30); res["C_measured_power_law_slopes_15_30gpc"]={f"m={m+1}":float(np.polyfit(np.log(r[sel]),np.log(A[sel,m]),1)[0]) for m in range(6)}
    json.dump(res,open("web_shells.json","w"),indent=1); json.dump(curves,open("cumulative_curves.json","w"))
    print(json.dumps({k:v for k,v in res.items() if k!="C_falloff_unit_k"},indent=1))
