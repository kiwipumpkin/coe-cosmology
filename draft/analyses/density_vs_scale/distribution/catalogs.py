# Number density of observed voids / clusters / superclusters vs radius, plus LCDM model curves.
import sys, json, numpy as np, pandas as pd
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
from scipy.integrate import quad
h=0.7; Om=0.31; rhom_h=Om*2.775e11          # h^2 Msun/Mpc^3
def rd(f):
    L=[l for l in open(f) if not l.startswith('#') and l.strip()]
    hdr=L[0].rstrip('\n').split('\t'); return pd.DataFrame([l.rstrip('\n').split('\t') for l in L[3:]],columns=hdr)
def dc(z): return 2997.92*quad(lambda x:1/np.sqrt(Om*(1+x)**3+1-Om),0,z)[0]      # Mpc/h
def shellvol(area_deg2,d1,d2): return area_deg2/41252.96*4/3*np.pi*(d2**3-d1**3)  # (Mpc/h)^3
G3=lambda V_h: V_h/h**3/1e9          # (Mpc/h)^3 -> Gpc^3
edges=np.logspace(0,3,31); cen=np.sqrt(edges[1:]*edges[:-1]); dlog=np.log10(edges[1]/edges[0])
def nd(R_Mpc,Vgpc3):
    c,_=np.histogram(R_Mpc,edges); return c/Vgpc3/dlog, c
out={}
# --- voids: Mao+2017 BOSS DR12 (cleaned catalog)
m=rd("data/mao2017_voids.tsv")
for c in ["z","Reff","delmin"]: m[c]=pd.to_numeric(m[c])
areas={"CMASS North":6851,"CMASS South":2525,"LOWZ North":5836,"LOWZ South":2501}
zr={"CMASS":(0.43,0.7),"LOWZ":(0.2,0.43)}
V=sum(shellvol(a,dc(zr[s.split()[0]][0]),dc(zr[s.split()[0]][1])) for s,a in areas.items())
n,c=nd(m.Reff/h,G3(V)); out["voids_Mao2017"]=dict(R=cen.tolist(),n_per_Gpc3_dex=n.tolist(),counts=c.tolist(),Vgpc3=G3(V),N=len(m),
    delmin_median=float(m.delmin.median()),delmin_q=[float(q) for q in m.delmin.quantile([0.1,0.9])])
# --- clusters: redMaPPer SDSS DR8, 0.1<z<0.3, lambda>=20; M200m from Simet+2017
r=rd("data/redmapper_dr8.tsv")
for c in ["zlambda","lambda"]: r[c]=pd.to_numeric(r[c],errors="coerce")
s=r[(r.zlambda>0.1)&(r.zlambda<0.3)]
M=10**14.344*(s["lambda"]/40)**1.33                       # Msun/h
R200m=(3*M/(4*np.pi*200*rhom_h))**(1/3)/h                  # Mpc
V=shellvol(10400,dc(0.1),dc(0.3))
n,c=nd(R200m,G3(V)); out["clusters_redMaPPer"]=dict(R=cen.tolist(),n_per_Gpc3_dex=n.tolist(),counts=c.tolist(),Vgpc3=G3(V),N=len(s),delta=200/1-1)
# --- superclusters: Liivamagi+2012 SDSS DR7 main, D=5 threshold; equivalent-sphere radius
l=rd("data/liivamagi2012_main50.tsv"); l["Vol"]=pd.to_numeric(l.Vol)
Req=(3*l.Vol/(4*np.pi))**(1/3)/h
V=shellvol(7646,55,565)
n,c=nd(Req,G3(V)); out["superclusters_Liivamagi2012"]=dict(R=cen.tolist(),n_per_Gpc3_dex=n.tolist(),counts=c.tolist(),Vgpc3=G3(V),N=len(l))
# --- LCDM models
def growth(z):
    a=1/(1+z); f=lambda x: 1/(x*np.sqrt(Om/x**3+1-Om))**3
    Dz=2.5*Om*np.sqrt(Om/a**3+1-Om)*quad(f,0,a)[0]; D0=2.5*Om*quad(f,0,1)[0]; return Dz/D0
Rm=np.logspace(-0.3,3,200)        # Mpc
# Tinker+2008 halo mass function, Delta=200m, z=0.2; plotted vs R200m
Dg=growth(0.2); Ms=np.logspace(12.5,16,120)                        # Msun/h
RL=(3*Ms/(4*np.pi*rhom_h))**(1/3)/h                                  # Lagrangian radius, Mpc
sig=np.array([sigma_TH(x/1000) for x in RL])*Dg
A,a,b,cc=0.186*(1.2)**-0.14,1.47*(1.2)**-0.06,2.57*(1.2)**-0.011,1.19   # z=0.2 scaling (Tinker08 eq 5-8)
fs=A*((sig/b)**-a+1)*np.exp(-cc/sig**2)
dlnsinv=-np.gradient(np.log(sig),np.log(Ms))
dndlnM=fs*rhom_h/Ms*dlnsinv*h**3*1e9                                 # per Gpc^3
R200=(3*Ms/(4*np.pi*200*rhom_h))**(1/3)/h
dndlogR=dndlnM*np.log(10)*3                                          # dlnM = 3 dlnR
out["model_Tinker08_z0.2"]=dict(R=R200.tolist(),n_per_Gpc3_dex=dndlogR.tolist())
# Jennings+2013 Vdn void size function (linear threshold -2.71, voids with delta=-0.8 inside r), z=0.5
Dg=growth(0.5); dv,dcol=-2.71,1.686; Dd=abs(dv)/(dcol+abs(dv))
rL=np.logspace(0,2.5,150); sg=np.array([sigma_TH(x/1000) for x in rL])*Dg
x=Dd*sg/abs(dv); j=np.arange(1,200)[:,None]
f=2*np.sum(np.exp(-(j*np.pi*x)**2/2)*j*np.pi*x**2*np.sin(j*np.pi*Dd),axis=0)
rr=rL*(1-0.8)**(-1/3); Vr=4/3*np.pi*rr**3
dlnsinv=-np.gradient(np.log(sg),np.log(rL))
dndlnr=f/Vr*dlnsinv*1e9
out["model_Jennings13_Vdn_z0.5"]=dict(R=rr.tolist(),n_per_Gpc3_dex=(dndlnr*np.log(10)).tolist())
json.dump(out,open("object_counts.json","w"),indent=1)
for k,v in out.items():
    if "N" in v: print(k,"N=",v["N"],"V=%.2f Gpc^3"%v["Vgpc3"], "total n=%.1f /Gpc^3"%(v["N"]/v["Vgpc3"]))
print("void delmin median",out["voids_Mao2017"]["delmin_median"],out["voids_Mao2017"]["delmin_q"])
