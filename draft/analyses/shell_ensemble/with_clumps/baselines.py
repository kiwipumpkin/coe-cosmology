"""Baselines with k=1..4 of the fitted four clumps (scan.json sym N=4, the 83x fit), with/without density-gradient compensation."""
import json, os, itertools, numpy as np
from scipy.optimize import nnls, least_squares
HERE=os.path.dirname(os.path.abspath(__file__))
g={"__file__":os.path.join(os.path.dirname(HERE),"surround","ensemble_surround.py"),"__name__":"imp"}
exec(open(g["__file__"]).read().split('\nif __name__')[0],g); os.chdir(HERE)
field3,amps,R,rms,RV,FLOOR,phi,KG,coe_vec,A_IN,S_MAX,L=[g[k] for k in ("field3","amps","R","rms","RV","FLOOR","phi","KG","coe_vec","A_IN","S_MAX","L")]
field=g["g"]["field"]
sc=json.load(open("../../scan.json"))["sym"]["4"]; Dc,k1,o1,k2,o2=sc["p"]
ang=np.radians([phi+o1,phi-o1,phi+o2,phi-o2]); kk=np.array([k1,k1,k2,k2])
P4=np.c_[Dc*np.cos(ang),Dc*np.sin(ang),np.zeros(4)]          # R_vis units, ring plane
F4=np.array([field3(a,0.0,Dc)[0] for a in ang])                 # unit-k patterns
C=coe_vec(phi+180.0); rC=np.linalg.norm(P4-C,axis=1)*RV
print("check 4-clump RMS/floor:",rms(R-kk@F4)/FLOOR,"(scan.json",sc["rms"]/FLOOR,")")

# ---- compensation: remove the clumps' mass from the surround region, spread ∝ rho(r_CoE) (so the profile shape is kept: rho -> (1-eps) rho)
rng=np.random.default_rng(12345); NS=400000
u=rng.normal(size=(NS*3,3)); u/=np.linalg.norm(u,axis=1)[:,None]; s=rng.uniform(1,S_MAX**3,NS*3)**(1/3); Q=s[:,None]*u
rq=np.linalg.norm(Q-C,axis=1); ok=(rq>=A_IN)&(rng.uniform(size=len(rq))<np.exp(-(rq-A_IN)/L)); Q=Q[ok][:NS]
Dq=np.linalg.norm(Q,axis=1); Fq=np.zeros_like(R)
for i in range(0,NS,20000):
    sl=slice(i,i+20000); Fq+=field3(np.arctan2(Q[sl,1],Q[sl,0]),np.arcsin(Q[sl,2]/Dq[sl]),Dq[sl]).sum(0)
Fcomp_unit=Fq/NS          # pattern per unit k of mass spread ∝ rho over the region
# region mass & eps
GPC=3.0857e25; RHO=2.9e-27
Vreg=(4/3*np.pi*(S_MAX**3-1))*RV**3*ok[:NS*3].mean()/np.exp(0)  # rough (acceptance includes gradient ~1)
def model(sub,comp,kvec=None):
    k=kk[list(sub)] if kvec is None else kvec
    w=k@F4[list(sub)] if len(sub) else 0*R
    if comp: w=w-np.sum(k)*Fcomp_unit
    return w
out=dict(clumps=[dict(idx=i,M_kg=float(kk[i]*KG),ell_deg=float(np.degrees(ang[i])%360),dist_from_us_gpc=float(Dc*RV),dist_from_CoE_gpc=float(rC[i]),
              zone=("inside D-d (CoE interior)" if rC[i]<43.9 else ("in shell" if rC[i]<=72.5 else "beyond D+d"))) for i in range(4)],
         base0=rms(R)/FLOOR,comp_unit_rms_over_floor_per_7p5e51=float(rms(Fcomp_unit)*7.5e51/KG/FLOOR),baselines={})
for k in range(1,5):
    rows=[]
    for sub in itertools.combinations(range(4),k):
        r_unc=rms(R-model(sub,False)); r_cmp=rms(R-model(sub,True))
        x,_=nnls(F4[list(sub)].T,R); r_refitM=rms(R-x@F4[list(sub)])
        rows.append(dict(sub=list(sub),rms_uncomp=r_unc/FLOOR,rms_comp=r_cmp/FLOOR,rms_massrefit=r_refitM/FLOOR,massrefit_kg=(x*KG).tolist()))
    best=min(rows,key=lambda d:d["rms_comp"])
    out["baselines"][k]=dict(best=best,all=rows)
    print(k,"best subset",best["sub"],f"uncomp {best['rms_uncomp']:.1f}  comp {best['rms_comp']:.1f}  (masses refit on same subset: {best['rms_massrefit']:.1f}; best any subset refit {min(r['rms_massrefit'] for r in rows):.1f})")
    for r in rows: print("    ",r["sub"],f"{r['rms_uncomp']:.1f} {r['rms_comp']:.1f} {r['rms_massrefit']:.1f}")
# free refit of k clumps (angles+masses, shared D, in ring plane), multistart, cheap
def gen(p,k):
    D=p[0]; return sum(p[1+2*j]*field(np.radians(p[2+2*j]),D) for j in range(k))
for k in (1,2,3):
    best=None
    for D0 in (1.2,1.5,1.7,2,3,5):
        for sd in range(12):
            r=np.random.default_rng(sd); p0=[D0]+sum([[r.normal(0,4e-3),r.uniform(0,360)] for _ in range(k)],[])
            f=least_squares(lambda p:gen(p,k)-R,p0,bounds=([1.05]+[-np.inf,-360]*k,[20]+[np.inf,720]*k),max_nfev=3000)
            if best is None or f.cost<best.cost: best=f
    p=best.x; out["baselines"][k]["free_refit"]=dict(rms=rms(R-gen(p,k))/FLOOR,D_gpc=p[0]*RV,M_kg=[p[1+2*j]*KG for j in range(k)],ell=[p[2+2*j]%360 for j in range(k)])
    print(k,"free refit",out["baselines"][k]["free_refit"])
np.savez("baseline_patterns.npz",F4=F4,kk=kk,Fcomp_unit=Fcomp_unit,P4=P4)
json.dump(out,open("baselines.json","w"),indent=1)
print(json.dumps(out["clumps"],indent=0)); print("comp pattern rms per 7.5e51 kg:",out["comp_unit_rms_over_floor_per_7p5e51"])
