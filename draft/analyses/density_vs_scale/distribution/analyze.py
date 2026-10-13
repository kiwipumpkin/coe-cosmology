import sys, json, numpy as np
from scipy.stats import norm
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
import lcdm
from lcdm import sigma_TH
from scipy.integrate import simpson
h=0.7
def sig_TH_G(R_Mpc,sg_h=4.0):   # linear sigma, top-hat R convolved with 4 Mpc/h Gaussian (2M++ pre-smoothing)
    Rh=R_Mpc*h; k=lcdm.k
    W=lcdm.W_th(k*Rh)*np.exp(-(k*sg_h)**2/2)
    return np.sqrt(lcdm.A*simpson(lcdm.P0*W**2*k**2/(2*np.pi**2),x=k))
# --- fine grid of radii from the 2M++ map
d=np.load("../data/twompp_density.npy").astype(np.float64)
n=257; dx=400/256; ax=(np.arange(n)-128)*dx
X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij"); r=np.sqrt(X**2+Y**2+Z**2)
F=np.fft.rfftn(d); bstar=0.3**0.55/0.43; RMAX=150.0
Rgrid=np.array([10,14,20,28,40,56,80,110,150.0])
lbins=np.linspace(-4,1.5,56)   # log10|delta| bins
rows=[]; H_pos=[]; H_neg=[]
for R in Rgrid:
    Rh=R*h; k=(r<=Rh).astype(float); k/=k.sum()
    a=(np.fft.irfftn(F*np.fft.rfftn(np.fft.ifftshift(k)),s=d.shape,axes=(0,1,2))/bstar)[r+Rh<=RMAX]
    w=np.ones_like(a)/a.size
    H_pos.append(np.histogram(np.log10(a[a>0]),lbins,weights=w[a>0])[0]); H_neg.append(np.histogram(np.log10(-a[a<0]),lbins,weights=w[a<0])[0])
    rows.append(dict(R=float(R),std=float(a.std()),mean=float(a.mean()),n_indep=float((RMAX/Rh)**3),
                     sig_lin=float(sigma_TH(R/1000)),sig_lin_presmoothed=float(sig_TH_G(R)),
                     q=[float(x) for x in np.percentile(a,[0.13,2.3,15.9,50,84.1,97.7,99.87])]))
    print(rows[-1]["R"],round(rows[-1]["std"],4),round(rows[-1]["sig_lin"],4),round(rows[-1]["sig_lin_presmoothed"],4),round(rows[-1]["n_indep"],1))
np.savez("heat_2mpp.npz",R=Rgrid,lbins=lbins,H_pos=np.array(H_pos),H_neg=np.array(H_neg))
# --- width trend: power-law fit of std vs R over 20-80 Mpc (exclude 10-14: pre-smoothing; >=110: <3 independent spheres)
sel=[r_ for r_ in rows if 20<=r_["R"]<=80]
x=np.log10([r_["R"] for r_ in sel]); y=np.log10([r_["std"] for r_ in sel])
X1=np.c_[np.ones(len(x)),x]; b=np.linalg.lstsq(X1,y,rcond=None)[0]; res=y-X1@b
cov=(res@res/(len(x)-2))*np.linalg.inv(X1.T@X1)
fitw=dict(logA=float(b[0]),beta=float(-b[1]),se_beta=float(np.sqrt(cov[1,1])),cov=cov.tolist(),points=[r_["R"] for r_ in sel])
print("width fit: std = 10^%.3f R^-%.3f ± %.3f"%(b[0],-b[1],np.sqrt(cov[1,1])))
# ratio data/LCDM(presmoothed) in the fit range -> effective bias-like factor
ratio=np.mean([r_["std"]/r_["sig_lin_presmoothed"] for r_ in sel]); print("2M++ std / LCDM linear (pre-smoothed) over 20-80 Mpc:",round(ratio,3))
# --- rarity of the required field
Rvis=14300.0
peak={1000:0.068,2000:0.055,4000:0.059}; rms={1000:0.0105,2000:0.0124,4000:0.0141}
rar={}
for R in [1000,2000,4000]:
    xx=np.log10(R); yy=b[0]+b[1]*xx; sy=np.sqrt(cov[0,0]+xx*xx*cov[1,1]+2*xx*cov[0,1])
    s_fit=10**yy; s_lo,s_hi=10**(yy-sy),10**(yy+sy); s_l=sigma_TH(R/1000)
    Nind=(Rvis/R)**3
    def stats(s):
        ns=peak[R]/s; p2=2*norm.sf(ns); return dict(sigma=float(s),n_sigma=float(ns),p_two_sided=float(p2),log10p=float((np.log(2)+norm.logsf(ns))/np.log(10)),
                                                  expected_in_visible_sphere=float(Nind*p2),rms_ratio=float(rms[R]/s))
    rar[R]=dict(N_indep_visible=float(Nind),trend=stats(s_fit),trend_wide=stats(s_hi),trend_narrow=stats(s_lo),lcdm=stats(s_l))
    for k_,v in rar[R].items():
        if isinstance(v,dict): print(R,k_,{kk:(f"{vv:.3g}") for kk,vv in v.items()})
    print(R,"N_indep",round(Nind,1))
json.dump(dict(rows=rows,width_fit=fitw,ratio_2mpp_to_lcdm=float(ratio),rarity=rar),open("distribution_stats.json","w"),indent=1)
