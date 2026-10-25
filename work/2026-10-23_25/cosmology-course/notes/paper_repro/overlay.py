import numpy as np,json,os
N='/workspace/cosmology-course/notes/'; OUT=N+'paper_repro/'
Q=json.load(open(N+'catwise-density-levelset-fit.json'))['migkas_fits']['soft_fixed15.2']
src=open(N+'catwise-density-levelset-migkas.py').read(); n2={'__name__':'x'}; exec(src[:src.index('def main')],n2)
ar=n2['_pc_args']('soft',Q['DR'],Q['LD'],Q['shape']); rq=n2['pc'].run('soft',ar,Q['mc'],0.0); af=rq['af']/rq['af'].mean()
l0=Q['l0']; bb,Tsp,_=np.load(OUT+'split_T.npy'); dT=(Tsp-Tsp.mean())/Tsp.mean()
A_beta=af  # af index = β (θ=0 away from CoE), so A(β)=af[β]
Y=n2['MA']; thd=(n2['LDEG']-l0)%360
r=np.corrcoef(A_beta,dT)[0,1]; pA=int(np.argmax(A_beta)); pT=int(bb[np.argmax(dT)])
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(16,10)); ax.errorbar(thd,Y,yerr=[n2['MLO'],n2['MHI']],fmt='o',color='k',ms=7,capsize=3,label='Migkas A/A_all (71 bins)')
ax.plot(np.arange(360),A_beta,lw=4,color='tab:orange',label='soft-core model A/A_all (R² 0.743)')
ax.set_xlim(0,360); ax.set_xlabel('β (deg, 0 = away from CoE)   [β = ℓ − 15.2°]',fontsize=24); ax.set_ylabel('A / A_all (dimensionless)',fontsize=24,color='tab:orange'); ax.tick_params(labelsize=20)
t=ax.twinx(); t.plot(bb,dT*1e5,lw=4,ls='--',color='tab:blue',label='split-shell CMB ΔT/T'); t.set_ylabel('CMB ΔT/T  (×10⁻⁵)',fontsize=24,color='tab:blue'); t.tick_params(labelsize=20)
sec=ax.secondary_xaxis('top',functions=(lambda x:x+l0,lambda x:x-l0)); sec.set_xlabel('Galactic ℓ (deg)',fontsize=22); sec.tick_params(labelsize=18)
h1,l1=ax.get_legend_handles_labels(); h2,l2=t.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,fontsize=18,loc='lower center',bbox_to_anchor=(0.5,0.02))
ax.set_title(f'A/A_all and split-shell CMB ΔT/T; correlation r = {r:.2f}',fontsize=24)
fig.tight_layout(); fig.savefig(OUT+'overlay_A_dT.png',dpi=100)
print(dict(r=r,A_peak_beta=pA,A_peak_l=(pA+l0)%360,dT_peak_beta=pT,dT_peak_l=(pT+l0)%360,A_min_beta=int(np.argmin(A_beta)),dT_amp_1e5=float(dT.max()*1e5)))
