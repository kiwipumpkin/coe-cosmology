import json,os,numpy as np
N='/workspace/cosmology-course/notes/'
src=open(N+'catwise-density-levelset-migkas.py').read(); ns={'__name__':'x'}; exec(src[:src.index('def main')],ns)
os.chdir(N+'paper_repro'); pc=ns['pc']; import bubble_core as bc
p=json.load(open(N+'catwise-density-levelset-fit.json'))['migkas_fits']['soft_fixed15.2']
a=ns['_pc_args']('soft',p['DR'],p['LD'],p['shape']); rr=pc.run('soft',a,p['mc'],0.0)
q,rho=bc.geometry(bc.TH1,p['DR']); X,V,_=pc.integrate(q,rr['mu'],pc.RG,rr['MT'],0.0,rr['Rb'],0.01,200)
DG=p['D_Gpc']; t0=rr['Ht']*bc.TH0; C=0.306601*t0/DG
af=np.hypot(*(V[1:]-V[0]).T)/rho; assert np.allclose(af,rr['af'])      # same velocities the fit uses
Q=q*DG; d=rho*DG; U=np.hypot(*V[0])/C; vin=np.hypot(*V[1+180])/C; vout=np.hypot(*V[1])/C
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(18,13)); th=np.linspace(0,2*np.pi,500)
g=np.linspace(-5,72,400); xx,yy=np.meshgrid(g,np.linspace(-24,24,300)); S=np.interp(np.hypot(xx,yy)/DG,pc.RG,rr['S'])
ax.pcolormesh(g,np.linspace(-24,24,300),np.log(S),cmap='Greys',alpha=.3,shading='auto')
ax.plot(0,0,'*',color='r',ms=40,zorder=6); ax.text(0,-3.5,f'CoE + mass\nM(<D) = {p["mc"]:.2f} M_crit\n(soft core, L = {p["L_Gpc_locked"]/1e3:.1f}×10³ Gpc)',ha='center',va='top',fontsize=20)
for r,c_,lab in ((DG-d,'tab:blue',f'shell D−d = {DG-d:.1f}'),(DG+d,'tab:red',f'shell D+d = {DG+d:.1f}')):
    ax.plot(r*np.cos(th),r*np.sin(th),color=c_,lw=3.5)
    for ang in (22,-22): 
        u=np.array([np.cos(np.radians(ang)),np.sin(np.radians(ang))]); ax.annotate('',xy=(r+4)*u,xytext=r*u,arrowprops=dict(arrowstyle='-|>',lw=3,color=c_))
    ax.text(r*np.cos(np.radians(-30))+1,r*np.sin(np.radians(-30)),lab+' Gpc',color=c_,fontsize=20,ha='center')
ax.plot(Q[1:,0],Q[1:,1],color='green',lw=4); ax.plot(*Q[0],'o',color='orange',ms=24,mec='k',zorder=7)
ax.text(Q[0,0],-3,f'us, D = {DG:.1f} Gpc',ha='center',va='top',fontsize=22)
SG=1.5   # Gpc per c (global); SP: Gpc per c (relative)
for P,v,c_,lab in ((Q[0],V[0],'orange','U'),(Q[1+180],V[1+180],'tab:blue','v_in'),(Q[1],V[1],'tab:red','v_out')):
    tip=P+SG*v/C; ax.annotate('',xy=tip,xytext=P,arrowprops=dict(arrowstyle='-|>',lw=5,color=c_),zorder=8)
    ax.text(tip[0],tip[1]+1.2,f'{lab} = {np.hypot(*v)/C:.1f}c',color=c_,fontsize=20,ha='center',fontweight='bold',bbox=dict(fc='w',ec='none',alpha=.8))
idx=np.arange(0,360,30); vrel=(V[1:]-V[0])[idx]/C; SP=3.0
for i,(P,w) in zip(idx,zip(Q[1:][idx],vrel)):
    ax.annotate('',xy=P+SP*w,xytext=P,arrowprops=dict(arrowstyle='-|>',lw=3,color='k'))
for i in (45,135,225,315):
    ax.annotate('',xy=Q[0]+(Q[1+i]-Q[0])*0.12,xytext=Q[1+i],arrowprops=dict(arrowstyle='->',lw=2.5,color='gold',ls='--'))
ax.text(Q[0,0]+d+1,d*0.8,f'visible circle d = {d:.1f} Gpc',color='green',fontsize=20)
ax.text(-4,22,f'Key: coloured arrows = speed from CoE, {SG} Gpc per c\nblack = velocity relative to us (sets A/A_all), {SP} Gpc per c ({np.abs(vrel).max():.1f}c max)\ngold dashed = light from circle to us;  θ=0 (right) is ℓ₀ = {p["l0"]:.1f}°',fontsize=18,va='top',bbox=dict(fc='w',ec='k'))
ax.set_xlim(-5,72); ax.set_ylim(-24,24); ax.set_aspect('equal')
ax.set_xlabel('Gpc from CoE',fontsize=22); ax.set_ylabel('Gpc',fontsize=22); ax.tick_params(labelsize=18)
ax.set_title('Soft-core best fit today (paper code geometry and velocities)',fontsize=26)
fig.tight_layout(); fig.savefig('schematic_fitted.png',dpi=90)
print(dict(D=DG,d=d,U_c=U,vin_c=vin,vout_c=vout,vrel_max_c=float(np.hypot(*((V[1:]-V[0])/C).T).max()),pos_today_int=np.hypot(*X[0])*DG))
