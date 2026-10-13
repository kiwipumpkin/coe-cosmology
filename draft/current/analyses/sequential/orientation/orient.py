"""Orientation of cluster relative to its void in the fitted near arrangements (deep voids, Scenario A)."""
import sys, os, json, numpy as np
OD=os.path.dirname(os.path.abspath(__file__)); SEQ=os.path.dirname(OD); sys.path.insert(0,SEQ)
from seq_core import *
DELTA_C=100.0
def clu_R(M,x): return float((3*abs(M)/(4*np.pi*DELTA_C*rho_gpc(np.asarray(x).reshape(1,3))[0]*GPC**3))**(1/3))
lf=json.load(open(os.path.join(SEQ,"local_fit.json")))
def unit(v): return v/np.linalg.norm(v)
def ringd(x): return float(np.sqrt((np.hypot(x[0],x[1])-RMIN)**2+x[2]**2))
rows=[]
for run in lf["runs"]:
    t=run["traj"]; vs=[p for p in t if p["kind"]=="void"]; cs=[p for p in t if p["kind"]=="cluster"]
    V=np.array([p["x"] for p in vs])
    for k,(v,c) in enumerate(zip(vs,cs)):
        xv=np.array(v["x"]); xc=np.array(c["x"]); o=xc-xv; L=np.linalg.norm(o); u=o/L
        rhat=unit(xv); chat=unit(C_G-xv); phat=unit(np.cross([0,0,1.0],xv))   # tangential in ring plane
        nearest=np.argmin(np.linalg.norm(V-xc,axis=1))
        rows.append(dict(seed=run["seed"],k=k,Rv=v["R"],Rc=c["R"],L=L,L_Rv=L/v["R"],
            cos_radial=float(u@rhat),            # +1 cluster directly behind its void (farther from us), -1 on our side
            cos_coe=float(u@chat),               # +1 cluster toward the CoE from its void
            cos_tan=float(abs(u@phat)),cos_z=float(abs(u[2])),
            d_void=float(np.linalg.norm(xv)),d_clu=float(np.linalg.norm(xc)),ring_void=ringd(xv),ring_clu=ringd(xc),
            void_angle_from_coe=float(np.degrees(np.arccos(unit(xv)@unit(C_G)))),clu_angle_from_coe=float(np.degrees(np.arccos(unit(xc)@unit(C_G)))),
            void_z=float(xv[2]),clu_z=float(xc[2]),az_sep=float(np.degrees(np.arccos(np.clip(unit(xv[:2])@unit(xc[:2]),-1,1)))),
            partner_is_nearest_void=bool(nearest==k),mass_ratio=float(c["M"]/-v["M"]),gain=None))
json.dump(rows,open(os.path.join(OD,"pairs.json"),"w"),indent=1)
a=lambda k:np.array([r[k] for r in rows])
S=dict(n=len(rows),
 offset_gpc_median=float(np.median(a("L"))),offset_gpc_p10_p90=[float(np.percentile(a("L"),10)),float(np.percentile(a("L"),90))],
 offset_in_Rv_median=float(np.median(a("L_Rv"))),
 frac_cluster_nearer_us=float(np.mean(a("d_clu")<a("d_void"))),frac_cos_radial_neg=float(np.mean(a("cos_radial")<0)),
 median_cos_radial=float(np.median(a("cos_radial"))),
 frac_cluster_toward_coe=float(np.mean(a("cos_coe")>0)),median_cos_coe=float(np.median(a("cos_coe"))),
 median_abs_cos_tangential=float(np.median(a("cos_tan"))),median_abs_cos_vertical=float(np.median(a("cos_z"))),
 frac_mostly_tangential=float(np.mean(a("cos_tan")>np.maximum(np.abs(a("cos_radial")),a("cos_z")))),
 frac_mostly_radial=float(np.mean(np.abs(a("cos_radial"))>np.maximum(a("cos_tan"),a("cos_z")))),
 frac_mostly_vertical=float(np.mean(a("cos_z")>np.maximum(a("cos_tan"),np.abs(a("cos_radial"))))),
 median_void_dist=float(np.median(a("d_void"))),median_clu_dist=float(np.median(a("d_clu"))),
 median_void_ring_dist=float(np.median(a("ring_void"))),median_clu_ring_dist=float(np.median(a("ring_clu"))),
 frac_clu_closer_to_ring=float(np.mean(a("ring_clu")<a("ring_void"))),
 median_abs_z_void=float(np.median(np.abs(a("void_z")))),median_abs_z_clu=float(np.median(np.abs(a("clu_z")))),
 median_azimuth_sep_deg=float(np.median(a("az_sep"))),
 void_angle_from_coe=np.percentile(a("void_angle_from_coe"),[10,50,90]).tolist(),clu_angle_from_coe=np.percentile(a("clu_angle_from_coe"),[10,50,90]).tolist(),
 frac_void_in_coe_hemisphere=float(np.mean(a("void_angle_from_coe")<90)),frac_clu_in_coe_hemisphere=float(np.mean(a("clu_angle_from_coe")<90)),
 frac_partner_is_nearest_void=float(np.mean(a("partner_is_nearest_void"))),mass_ratio=float(np.median(a("mass_ratio"))))
print(json.dumps(S,indent=1)); json.dump(S,open(os.path.join(OD,"stats.json"),"w"),indent=1)
