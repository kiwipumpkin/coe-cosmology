import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
D=json.load(open('l0scan_results.json'));R0=json.load(open('interp_results.json'))
T={}
for d in D:T.setdefault(d['run'],{})[d['l0']]=d
def get(k):
    l=sorted(T[k]);return np.array(l),np.array([T[k][x]['chi2'] for x in l]),np.array([T[k][x]['R2'] for x in l])
def panel(keys,fn,title):
    fig,ax=plt.subplots(1,2,figsize=(13,5));cm=plt.cm.viridis(np.linspace(0,1,len(keys)))
    for c,k in zip(cm,keys):
        l,x,r=get(k);ax[0].plot(l,x,'o-',ms=3,color=c,label=k);ax[1].plot(l,r,'o-',ms=3,color=c,label=k)
    ax[0].set_ylabel('chi2 (71 bins)');ax[1].set_ylabel('R2')
    for a in ax:a.set_xlabel('l0 fixed (deg)');a.grid(alpha=.3)
    ax[0].legend(fontsize=7);fig.suptitle(title);fig.tight_layout();fig.savefig(fn,dpi=90)
path=sorted(k for k in T if k.startswith('path'));hyb=sorted(k for k in T if k.startswith('hyb'))
panel(path,'l0scan_path.png','Migkas fit with l0 fixed: path steps 0 (REF) to 10 (RED)')
panel(hyb,'l0scan_hybrid.png','Migkas fit with l0 fixed: hybrids')
C3=[0,0.02,0.04,0.06,0.08,0.1,0.12,0.1353];FR=[30,100,300,1000,3000,10000,20000]
M=np.full((len(C3),len(FR)),np.nan);A=np.full_like(M,np.nan)
for i,c in enumerate(C3):
    for j,f in enumerate(FR):
        k='grid_C3_%g_f%d'%(c,f)
        if k in T:l,x,r=get(k);M[i,j]=x.min();A[i,j]=l[x.argmin()]
fig,ax=plt.subplots(figsize=(9,6));im=ax.imshow(M,origin='lower',cmap='viridis',aspect='auto');plt.colorbar(im,label='min chi2 over l0 in [-10,100]')
for i in range(len(C3)):
    for j in range(len(FR)):
        if np.isfinite(M[i,j]):ax.text(j,i,'%.2f\nl0=%g'%(M[i,j],A[i,j]),ha='center',va='center',fontsize=7,color='w')
ax.set_xticks(range(len(FR)),FR);ax.set_yticks(range(len(C3)),C3);ax.set_xlabel('frac');ax.set_ylabel('C3');fig.tight_layout();fig.savefig('l0scan_grid.png',dpi=90)
print('| run | ratio | '+' | '.join('l0=%d'%x for x in [0,4,6,10,15,20])+' | min chi2 | at l0 | R2 at min |')
for k in path+hyb:
    l,x,r=get(k);i=x.argmin()
    print('| %s | %.5f | '%(k,R0[k]['ratio'])+' | '.join('%.2f'%T[k][v]['chi2'] if v in T[k] else '-' for v in [0,4,6,10,15,20])+' | %.2f | %g | %.3f |'%(x[i],l[i],r[i]))
print('grid min chi2 / argmin l0 (rows C3, cols frac):')
for i,c in enumerate(C3):print(c,' '.join('%.2f@%g'%(M[i,j],A[i,j]) for j in range(len(FR))))
