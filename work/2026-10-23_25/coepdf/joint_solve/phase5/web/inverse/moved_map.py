# (c) Which bodies moved (seed 1): clusters (delta>0, red) vs voids (delta<0, blue), from where to where, marker size ~ |M|.
# Panels: sky (Mollweide lon/lat of direction) original (open) -> final (filled); x-z projection with arrows; r_orig vs r_final; displacement vs |M|.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
W='/workspace/coepdf/joint_solve/phase5/web/inverse/';R=21.7439;z=np.load(W+'inv_s1.npz')
for k in ('y0','y','M'):
    if not np.all(np.isfinite(z[k])):raise SystemExit('NONFINITE '+k)
m=z['moved'];y0=z['y0'][m];y=z['y'][m];M=z['M'][m];cl=M>0;s=4+60*np.abs(M)/np.abs(z['M']).max()
def ang(v):r=np.linalg.norm(v,axis=1);return np.arctan2(v[:,1],v[:,0]),np.arcsin(v[:,2]/r)
fig=plt.figure(figsize=(15,11))
a=fig.add_subplot(221,projection='mollweide')
for sel,c,n in ((cl,'r','clusters'),(~cl,'b','voids')):
    l0,b0=ang(y0[sel]);l1,b1=ang(y[sel]);a.scatter(l0,b0,s=s[sel],facecolors='none',edgecolors=c,lw=0.4,alpha=0.5);a.scatter(l1,b1,s=s[sel],c=c,alpha=0.6,label=f'{n} final (n={sel.sum()})')
a.set_title('directions: open = original, filled = final (z axis = CoE axis frame of the run)');a.legend(fontsize=7,loc='lower right');a.grid(alpha=0.3)
a=fig.add_subplot(222)
for sel,c in ((cl,'r'),(~cl,'b')):
    for p,q in zip(y0[sel],y[sel]):a.annotate('',q[[0,2]],p[[0,2]],arrowprops=dict(arrowstyle='->',color=c,lw=0.4,alpha=0.5))
    a.scatter(y[sel][:,0],y[sel][:,2],s=s[sel],c=c,alpha=0.6)
t=np.linspace(0,2*np.pi,300);a.plot(R*np.cos(t),R*np.sin(t),'k',lw=0.8);a.plot(0.975*R*np.cos(t),0.975*R*np.sin(t),'k:',lw=0.6)
a.set_aspect('equal');a.set_xlabel('x [Gpc]');a.set_ylabel('z [Gpc]');a.set_title('x–z projection, arrows original → final; circle = CMB shell')
a=fig.add_subplot(223);r0=np.linalg.norm(y0,axis=1)/R;r1=np.linalg.norm(y,axis=1)/R
for sel,c,n in ((cl,'r','clusters'),(~cl,'b','voids')):a.scatter(r0[sel],r1[sel],s=s[sel],c=c,alpha=0.5,label=n)
a.plot([0,1],[0,1],'k:');a.axhline(0.975,color='gray',lw=0.6);a.set_xlabel('original r/R');a.set_ylabel('final r/R');a.legend(fontsize=8);a.set_title('radial origin vs destination')
a=fig.add_subplot(224);d=np.linalg.norm(y-y0,axis=1)
for sel,c,n in ((cl,'r','clusters'),(~cl,'b','voids')):a.scatter(np.abs(M[sel]),d[sel],s=6,c=c,alpha=0.5,label=n)
a.set_xscale('log');a.set_xlabel('|M| [Msun]');a.set_ylabel('displacement [Gpc]');a.legend(fontsize=8);a.set_title('displacement vs body mass')
fig.suptitle(f'Inverse problem seed 1: {m.sum()} of {len(m)} bodies moved ({cl.sum()} clusters, {(~cl).sum()} voids)');fig.tight_layout();fig.savefig(W+'moved_map.png',dpi=120)
print(dict(n=int(m.sum()),clusters=int(cl.sum()),voids=int((~cl).sum()),Mabs_med=float(np.median(np.abs(M))),disp_med_cl=float(np.median(d[cl])),disp_med_void=float(np.median(d[~cl])),r1_mean_cl=float(r1[cl].mean()),r1_mean_void=float(r1[~cl].mean())))
