# Copy of moved_map.py writing moved_map_s<seed>.png so the seed-1 moved_map.png is not overwritten.
# v3 (seed given, default 1): which bodies changed. Clusters (final delta>0, red) vs voids (blue); open = original, filled = final; size ~ |delta|.
# Panels: sky directions; x-z projection with arrows; r_orig vs r_final; |delta| start vs final. Usage: python3 moved_map.py [seed] [dir/]. Fails on NaN.
import sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
sd=sys.argv[1] if len(sys.argv)>1 else '1';W='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/'+(sys.argv[2] if len(sys.argv)>2 else '');R=21.7439;z=np.load(W+f'inv3_s{sd}.npz')
for k in ('y0','y','delta','delta0'):
    if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
m=z['changed'];y0=z['y0'][m];y=z['y'][m];d0=z['delta0'][m];d1=z['delta'][m];cl=d1>0;s=4+60*np.abs(d1)/max(np.abs(d1).max(),1e-30)
def ang(v):r=np.linalg.norm(v,axis=1);return np.arctan2(v[:,1],v[:,0]),np.arcsin(v[:,2]/r)
fig=plt.figure(figsize=(15,11));a=fig.add_subplot(221,projection='mollweide')
for q,c,nm in ((cl,'r','clusters'),(~cl,'b','voids')):
    l0,b0=ang(y0[q]);l1,b1=ang(y[q]);a.scatter(l0,b0,s=8,facecolors='none',edgecolors=c,lw=0.4,alpha=0.4);a.scatter(l1,b1,s=s[q],c=c,alpha=0.6,label=f'{nm} final (n={q.sum()})')
a.set_title('directions: open = original, filled = final (run frame)');a.legend(fontsize=7,loc='lower right');a.grid(alpha=0.3)
a=fig.add_subplot(222)
for q,c in ((cl,'r'),(~cl,'b')):
    for p,t in zip(y0[q],y[q]):a.annotate('',t[[0,2]],p[[0,2]],arrowprops=dict(arrowstyle='->',color=c,lw=0.4,alpha=0.5))
    a.scatter(y[q][:,0],y[q][:,2],s=s[q],c=c,alpha=0.6)
t=np.linspace(0,2*np.pi,300);a.plot(R*np.cos(t),R*np.sin(t),'k',lw=0.8);a.plot(0.975*R*np.cos(t),0.975*R*np.sin(t),'k:',lw=0.6)
a.set_aspect('equal');a.set_xlabel('x [Gpc]');a.set_ylabel('z [Gpc]');a.set_title('x-z projection, original -> final; solid = CMB shell, dotted = 0.975 R')
a=fig.add_subplot(223);r0=np.linalg.norm(y0,axis=1)/R;r1=np.linalg.norm(y,axis=1)/R
for q,c,nm in ((cl,'r','clusters'),(~cl,'b','voids')):a.scatter(r0[q],r1[q],s=s[q],c=c,alpha=0.5,label=nm)
a.plot([0,1],[0,1],'k:');a.axhline(0.975,color='gray',lw=0.6);a.set_xlabel('original r/R');a.set_ylabel('final r/R');a.legend(fontsize=8);a.set_title('radial origin vs destination')
a=fig.add_subplot(224)
for q,c,nm in ((cl,'r','final clusters'),(~cl,'b','final voids')):a.scatter(d0[q],d1[q],s=6,c=c,alpha=0.5,label=nm)
a.axhline(-1,color='k',ls=':',lw=0.8,label='delta = -1 (empty void)');a.set_yscale('symlog',linthresh=0.1);a.set_xlabel('start delta');a.set_ylabel('final delta');a.legend(fontsize=8);a.set_title('density contrast start vs final')
fig.suptitle(f'v3 inverse web seed {sd}: {m.sum()} of {len(m)} bodies changed ({z["moved"].sum()} moved, {(z["mch"]&~z["moved"]).sum()} mass-only)');fig.tight_layout();fig.savefig(W+f'moved_map_s{sd}.png',dpi=120)
