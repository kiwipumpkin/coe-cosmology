"""Self-consistent redo using Chris's spreadsheet definition (cosmology/papers/spreadsheet_model/columns.tex):
V = k r about the CoE (Sedov-Taylor homologous), rho ∝ exp(-r/L) from r=0 with M(<D) fixed, relative accel = a(P)-a(obs),
'recession v with gravity' = |v_rel + a_rel dt|. For V = k r, v_rel = k s n (radial from us), so the transverse term is zero and
d|v_rel|/dt = sign(k) a_los;  d2s/dt2 = a_los  (identical for k>0 and k<0).  Homologous V=kr is exactly consistent with gravity
only for uniform density (E(r) ∝ r^2); other L are the spreadsheet's own (instantaneous) combination."""
import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
D=58.2
th=np.radians(np.linspace(0,180,3601)); ct,st=np.cos(th),np.sin(th)
def Menc(r,L):  # exp sphere, normalised M(D)=1 ; L=0 -> point mass
    if L==0: return np.ones_like(r)
    m=lambda x:1-np.exp(-x)*(1+x+x*x/2); return m(r/L)/m(D/L)
def alos(s,L):
    rg=np.sqrt(D*D+s*s+2*D*s*ct); f=lambda r:Menc(r,L)/r**3
    return -f(rg)*(s+D*ct)+f(np.array(D))*D*ct   # units GM(D)/Gpc^2 per Gpc
def frac(m): return float(np.trapezoid(m*st,th)/2)
def band(y):
    i=np.where(np.diff(np.sign(y))!=0)[0]; return sorted(round(float(180-np.degrees(th[k])),1) for k in i)
out={}
for L in [0,2,4,6.07,10,20,50,922.65]:
  for s in [1,7,14]:
    y=alos(s,L); out[f"L={L}|s={s}"]=dict(frac_pos=round(frac(y>0),3),flip_deg_from_CoE=band(y))
json.dump(out,open("selfconsistent.json","w"),indent=1)
for k,v in out.items(): print(k,v)
lon=np.radians(np.linspace(-180,180,721)); lat=np.radians(np.linspace(-90,90,361)); LO,LA=np.meshgrid(lon,lat)
ang=np.degrees(np.arccos(np.clip(np.cos(LA)*np.cos(LO),-1,1)))
def mp(y,title,fn,lab):
    Y=np.interp(ang,180-np.degrees(th[::-1]),y[::-1]); m=np.abs(Y).max()
    fig=plt.figure(figsize=(9,5.2)); ax=fig.add_subplot(111,projection="mollweide")
    pc=ax.pcolormesh(LO,LA,Y,cmap="RdBu_r",norm=TwoSlopeNorm(0,-m,m),shading="auto")
    ax.contour(LO,LA,Y,levels=[0],colors="k",linewidths=1); ax.plot(0,0,"k*",ms=14); ax.text(0.06,0.06,"CoE",fontsize=11)
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.grid(alpha=.3); ax.set_title(title,fontsize=12,pad=12)
    cb=fig.colorbar(pc,ax=ax,orientation="horizontal",fraction=.05,pad=.06); cb.set_label(lab,fontsize=10)
    fig.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(fig)
y=alos(7,6.07)
mp(y,"Expanding (V = k r): relative acceleration along sightline\nexp profile L = 6.07 Gpc, s = 7 Gpc   (red = away from us)","sc_expanding_L6.07_s7.png","a_rel · n  (relative units)")
mp(y,"Collapsing (V = −k r): relative acceleration along sightline\nexp profile L = 6.07 Gpc, s = 7 Gpc   (red = away from us)","sc_collapsing_L6.07_s7.png","a_rel · n  (relative units)")
# claim (1) check: transverse term needed at th=90 for non-homologous v(r): need (v(rg) D/rg - v(D))^2/s > f(rg) s
