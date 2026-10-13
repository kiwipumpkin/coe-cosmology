import sys, json, numpy as np
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
h=0.7; rho_h=0.3*2.775e11      # mean matter density, h^2 Msun/Mpc^3 (Om=0.3)
def dm(M_h,V_h3): return M_h/(rho_h*V_h3)-1          # derived delta from mass [Msun/h] & volume [(Mpc/h)^3]
def Req(V_h3): return (3*V_h3/(4*np.pi))**(1/3)/h     # Mpc
TM=json.load(open("twompp_extremes.json"))
# name: (R_Mpc, |delta|, group, label?)  -- groups: def, mass, gal, tpp, void_gal, void_mass, disputed, tpp_min
P=[]
add=lambda *a: P.append(dict(zip(["name","R","d","group","sign","label"],a)))
add("Coma R200",2.84,645,"def",+1,"Coma")
add("Perseus R200",1.79,645,"def",+1,None)
add("Virgo zero-velocity",7.2,8.0e14/(0.3*2.775e11*h*h*4/3*np.pi*7.2**3)-1,"mass",+1,"Virgo")
add("Shapley core",11.4,20,"mass",+1,"Shapley core")
cr=[dm(0.6e16,4/3*np.pi*12.5**3),dm(12e16,4/3*np.pi*12.5**3)]
add("Corona Borealis SC",12.5/h,np.sqrt(cr[0]*cr[1]),"mass",+1,None)
add("Saraswati",36.6,1.62,"mass",+1,"Saraswati")
add("Shapley whole",Req(1.2e5),dm(1.8e16,1.2e5),"mass",+1,None)
add("Sloan Great Wall",Req(7.2e5),dm(1.2e17,7.2e5),"mass",+1,"Sloan Great Wall")
for n,M,V in [("A",2.1e16,67500),("B",1.1e16,70848),("C",0.7e16,19008),("D",0.3e16,13635)]:
    add(f"BOSS Great Wall {n}",Req(V),dm(M,V),"mass",+1,"BOSS Great Wall" if n=="A" else None)
add("Horologium-Reticulum SC",35,1.4,"gal",+1,None)
add("Vela SC",50,0.635,"gal",+1,None)
for r in TM["rows"]:
    if r["R_hMpc"]>=10:
        add(f"2M++ max R={r['R_hMpc']}h-1",r["R_Mpc"],r["max_dm"],"tpp",+1,None)
        add(f"2M++ min R={r['R_hMpc']}h-1",r["R_Mpc"],-r["min_dm"],"tpp_min",-1,None)
add("Local Void",23,0.8,"void_gal",-1,"Local Void")
add("Bootes void",44,0.75,"void_gal",-1,"Boötes")
add("BOSS largest void (central min)",64.05/h,0.89,"void_gal",-1,None)
add("SDSS voids lensing",40,0.4,"void_mass",-1,None)
add("Eridanus (Szapudi15)",220/h,0.14,"void_mass",-1,"Eridanus")
add("Eridanus (Kovacs&GB16, central)",195/h,0.25,"void_mass",-1,None)
add("Local Hole: Boehringer20 (100)",100,0.30,"disputed",-1,None)
add("Local Hole: Boehringer20 (140)",140,0.20,"disputed",-1,None)
add("Local Hole: Whitbourn&Shanks14",150/h,0.15,"disputed",-1,None)
add("Local Hole: Keenan13",300,0.33,"disputed",-1,"Local Hole claims")
add("Local Hole: Wong22",0.075*2998/h,0.21,"disputed",-1,None)
json.dump(P,open("points_v2.json","w"),indent=1,default=float)
def fit(sel,tf=lambda d:d):
    R=np.array([p["R"] for p in sel]); d=np.array([p["d"] for p in sel])
    x=np.log10(R); y=np.log10(tf(d)); n=len(x); X=np.c_[np.ones(n),x]
    b=np.linalg.lstsq(X,y,rcond=None)[0]; r=y-X@b; s2=r@r/(n-2); cov=s2*np.linalg.inv(X.T@X)
    return dict(logA=b[0],beta=-b[1],se_beta=float(np.sqrt(cov[1,1])),scatter_dex=float(np.sqrt(r@r/n)),n=n,cov=cov.tolist(),names=[p["name"] for p in sel])
G=lambda *gs:[p for p in P if p["group"] in gs]
F={}
F["over_all"]=fit(G("def","mass","gal","tpp"))
F["over_noDef"]=fit(G("mass","gal","tpp"))
F["over_lit_noDef"]=fit(G("mass","gal"))
F["over_lit_withDef"]=fit(G("def","mass","gal"))
F["over_2mpp"]=fit(G("tpp"))
F["under_all"]=fit(G("void_gal","void_mass","disputed","tpp_min"))
F["under_noDisp"]=fit(G("void_gal","void_mass","tpp_min"))
F["under_lit_noDisp"]=fit(G("void_gal","void_mass"))
F["under_lit_withDisp"]=fit(G("void_gal","void_mass","disputed"))
F["under_2mpp"]=fit(G("tpp_min"))
sat=lambda d:d/(1-d)
F["under_noDisp_sat"]=fit(G("void_gal","void_mass","tpp_min"),sat)
F["under_all_sat"]=fit(G("void_gal","void_mass","disputed","tpp_min"),sat)
lc={R:float(-np.log(sigma_TH(R*1.05/1000)/sigma_TH(R/1000))/np.log(1.05)) for R in [10,30,100,300,1000,2000,4000]}
def pred(f,R,satur=False):
    x=np.log10(R); y=f["logA"]-f["beta"]*x; c=np.array(f["cov"]); sy=np.sqrt(c[0,0]+x*x*c[1,1]+2*x*c[0,1])
    v,lo,hi=10**y,10**(y-sy),10**(y+sy)
    if satur: v,lo,hi=[q/(1+q) for q in (v,lo,hi)]
    return v,lo,hi
pk={1000:0.068,2000:0.055,4000:0.059}; rm={1000:0.0105,2000:0.0124,4000:0.0141}
for k,f in F.items():
    s=k.endswith("_sat")
    f["extrap"]={R:dict(zip(["fit","lo","hi"],pred(f,R,s))) for R in [830,1000,2000,3750,4000]}
    line=f"{k:20s} n={f['n']:2d} beta={f['beta']:.2f}±{f['se_beta']:.2f} scatter={f['scatter_dex']:.2f}dex |"
    for R in [1000,2000,4000]:
        v=f["extrap"][R]["fit"]; line+=f" {R}:{v:.1e} pk/fit={pk[R]/v:.3g} rms/fit={rm[R]/v:.2g} |"
    print(line)
print("LCDM beta (dln sigma/dln R):",{k:round(v,2) for k,v in lc.items()}, " white noise beta=1.5")
F["lcdm_beta"]=lc
json.dump(F,open("fits_v2.json","w"),indent=1,default=float)
for p in P: print(f"{p['name']:36s} R={p['R']:7.1f}  delta={p['sign']*p['d']:+.3g}  [{p['group']}]")
