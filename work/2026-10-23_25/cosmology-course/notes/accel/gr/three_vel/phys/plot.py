import pickle,json,numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.optimize import least_squares
exec(open('phys.py').read().split('fs=np.linspace')[0])
g=pickle.load(open('grid.pkl','rb'));fs,us,RT,CH=g['fs'],g['us'],g['RT'],g['CH']
i,j=np.unravel_index(np.nanargmin(CH),CH.shape);rat,u=RT[i,j],us[j];p=g['P'][(i,j)]
def st(m):return float(np.sum(((m-Y)/S)**2)),float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
cur={}
cur['best']=pred(rat,u,*p)
r=least_squares(lambda q:(pred(rat,u,q[0],q[1],264.021,48.253)-Y)/S,p[:2]);cur['planck']=pred(rat,u,r.x[0],r.x[1],264.021,48.253);pk=r.x
r=least_squares(lambda q:(pred(rat,u,q[0],0,264.021,48.253)-Y)/S,p[:1]);cur['f0']=pred(rat,u,r.x[0],0,264.021,48.253)
res={k:st(v) for k,v in cur.items()};res['planck_l0_bf']=list(pk);json.dump(res,open('curves.json','w'),indent=1,default=float);print(res)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':16})
fig,ax=plt.subplots(figsize=(13,8));ax.errorbar(L,Y,yerr=[MLO,MHI],fmt='o',color='k',ms=4,label='Migkas A/A_all (71 bins)')
lab={'best':'best fit, own-motion direction free','planck':'best fit, own motion toward Planck apex','f0':'frame speed 0, own motion toward Planck apex'}
o=np.argsort(L)
for k,c in zip(lab,['C0','C2','C3']):ax.plot(L[o],cur[k][o],c,lw=3,label=f"{lab[k]}: χ²={res[k][0]:.1f}, R²={res[k][1]:.2f}")
ax.set_xlabel('Galactic longitude (deg)');ax.set_ylabel('A / A_all');ax.legend(loc='lower center',bbox_to_anchor=(.5,1),frameon=False);fig.tight_layout();fig.savefig('phys_fit_migkas.png',dpi=120)
fig,ax=plt.subplots(figsize=(11,8));X,Yg=np.meshgrid(us,fs);d=CH-np.nanmin(CH)
cs=ax.contourf(X,Yg,np.minimum(d,20),levels=np.arange(0,21,1),cmap='viridis_r');plt.colorbar(cs,label='Δχ²')
ax.contour(X,Yg,d,levels=[2.3,6.2],colors='w');ax.set_xscale('log');ax.set_xlabel('expansion speed of space at our distance (units of c)');ax.set_ylabel('core mass / critical mass');ax.set_title('Δχ² (white: 1σ, 2σ). D is exactly degenerate (not shown)');fig.tight_layout();fig.savefig('degeneracy_core_vs_speed.png',dpi=120)
