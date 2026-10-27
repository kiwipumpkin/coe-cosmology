# Phase5 routes (NEW): present-day (t0) density and density at recombination vs r for the best fit of each route (+ previous unconstrained K-1 for reference).
# Reads fit_*.json; densities from homog.drec (same LTB physics as bg_mig/famP_mig). Output density_routes.png
import json,glob,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import homog as Hm
O=os.path.dirname(os.path.abspath(__file__))+'/'
xs=np.linspace(0.2,2.2,81);fig,ax=plt.subplots(1,2,figsize=(12,4.8))
runs=[('prev. unconstrained (K-1)',dict(D=27.737393,LD=0.42296181187611975,eta=3.068065))]+[(os.path.basename(f)[4:-5],json.load(open(f))) for f in sorted(glob.glob(O+'fit_[AB]*.json')) if 'ckpt' not in f]
for lab,r in runs:
    s=np.exp(-(xs/r['LD'])**r['eta']+(1/r['LD'])**r['eta']);dr=Hm.drec(r['D'],r['LD'],r['eta'],xs)[0]
    ax[0].semilogy(xs*r['D'],s,label=f"{lab}: D={r['D']:.2f}, eta={r['eta']:.3g}, L={r['LD']*r['D']:.3g}")
    ax[1].semilogy(xs*r['D'],np.where(dr>0,dr,np.nan),label=lab)
for a,t in zip(ax,['today: rho(t0,r)/rho(t0,D)','recombination: rho(t_rec,r)/rho(t_rec,D)']):
    a.set_xlabel('r from CoE [Gpc]');a.set_title(t);a.grid(alpha=.3);a.legend(fontsize=7)
plt.tight_layout();plt.savefig(O+'density_routes.png',dpi=130)
