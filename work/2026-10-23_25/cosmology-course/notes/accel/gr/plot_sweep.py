import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
S=pickle.load(open('sweep_static.pkl','rb')); M=pickle.load(open('sweep_mcv.pkl','rb'))
rows=[]
for v in sorted(S): to,res=S[v]; rows.append((v,'static',to,res[(0.1,0.9999)],res[(0.1,0)]))
for v in sorted(M): dt,res=M[v]; rows.append((v,'McVittie',dt,res[(0.1,0.9999)],res[(0.1,0.0)]))
f,ax=plt.subplots(figsize=(14,8.5)); V=np.array([r[0] for r in rows])
E=np.array([stats(norm(r[3]))[:2] for r in rows]); N=np.array([stats(norm(r[4]))[:2] for r in rows])
for arr,st,lab in ((E,'-','shell expanding ≈c'),(N,'--','shell not expanding / carried by flow')):
    lo=V<1; hi=V>=1
    for j,(c,nm) in enumerate((('C0','dipole'),('C3','quadrupole'))):
        ax.plot(V[lo],arr[lo,j],st,marker='o',color=c,lw=3,ms=8,label=f"{nm}, {lab}")
        ax.plot(V[hi],arr[hi,j],st,marker='s',color=c,lw=3,ms=8)
ax.axvline(1,color='k',lw=2); ax.text(0.93,ax.get_ylim()[1]*0.97,'through space',ha='right',va='top',fontsize=15); ax.text(1.07,ax.get_ylim()[1]*0.97,'carried by space',ha='left',va='top',fontsize=15)
ax.set_xscale('log'); ax.set_xlim(0.008,12); ax.set_ylim(bottom=0); ax.grid(alpha=.3,which='both'); ax.tick_params(labelsize=14)
ax.set_xticks([0.01,0.05,0.1,0.2,0.3,0.5,1,2,3,5,10]); ax.set_xticklabels(['0.01','0.05','0.1','0.2','0.3','0.5','1','2','3','5','10'])
ax.set_xlabel("speed (units of c), log scale. Left: start speed, static mass (exact GR). Right: recession speed, McVittie with H R/c = v/10.\nThe two halves are different models: the jump at c is not continuity. View time differs per point (earliest full sky on the left; 2R/c on the right).",fontsize=12)
ax.set_ylabel("size of pattern in ΔT/T",fontsize=16); ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=2,fontsize=13)
plt.tight_layout(); plt.savefig('speed_sweep.png',dpi=120)
for r,e,n in zip(rows,E,N): print(r[0],r[1],'view %.1f'%r[2],'exp: %.4f %.4f %.4f'%(*e,stats(norm(r[3]))[3]),'notexp: %.4f %.4f %.4f'%(*n,stats(norm(r[4]))[3]))
