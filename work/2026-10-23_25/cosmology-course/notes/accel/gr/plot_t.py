import sys,pickle; sys.path.insert(0,'.'); from grcore import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=pickle.load(open('speeds_t_2.3_2.4.pkl','rb')); to=2.3; x=np.degrees(th)
def get(mu,u,v0,m): w,g=R[(to,mu,u,v0)]; return w if m==0 else norm(g[0])
for c in sorted(k for k in R if k[0]==to):
    print(c,'weak %.4f %.4f %.4f'%tuple(np.array(stats(get(*c[1:],0)))[[0,1,3]]),'| GR %.4f %.4f %.4f'%tuple(np.array(stats(get(*c[1:],1)))[[0,1,3]]))
def plot(models,fn,title):
    f,ax=plt.subplots(figsize=(14,9 if len(models)>1 else 8.5))
    for m,lw,al,tag in models:
        for v0,col in ((0.1,'C0'),(0.3,'C3'),(0.5,'C2')):
            for u,st,l in ((0,'-','not expanding'),(0.9999,'--','expanding ≈c')):
                r=get(0.1,u,v0,m); d,q,_,_=stats(r)
                ax.plot(x,r,st,color=col,lw=lw,alpha=al,label=f"{tag}start {v0}c, {l}: dipole {d:.4f}, quad {q:.4f}")
    ax.axhline(0,color='gray'); ax.set_xlim(0,360); ax.grid(alpha=.3); ax.set_ylabel("ΔT/T",fontsize=16); ax.tick_params(labelsize=14)
    ax.set_xlabel(title,fontsize=15); ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=2,fontsize=11 if len(models)>1 else 12); plt.tight_layout(); plt.savefig(fn,dpi=120); plt.close()
plot([(0,3,1,'')],'speeds_weak_t2.3.png',"direction (deg); weak field; motion toward 0°, GM/(Rc²) = 0.1 at 180°, viewed at t = 2.3")
plot([(1,3,1,'')],'speeds_gr_t2.3.png',"direction (deg); exact GR; motion toward 0°, GM/(Rc²) = 0.1 at 180°, viewed at t = 2.3")
plot([(0,4,0.35,'weak, '),(1,2.5,1,'GR, ')],'speeds_weak_vs_gr_t2.3.png',"direction (deg); weak field = thick faded, exact GR = thin; GM/(Rc²) = 0.1, t = 2.3")
