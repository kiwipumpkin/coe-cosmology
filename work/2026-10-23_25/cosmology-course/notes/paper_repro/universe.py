import json,os,numpy as np
N='/workspace/cosmology-course/notes/'
src=open(N+'catwise-density-levelset-migkas.py').read(); ns={'__name__':'x'}; exec(src[:src.index('def main')],ns)
os.chdir('/workspace/cosmology-course/notes/paper_repro')
pc=ns['pc']; import bubble_core as bc
F=json.load(open(N+'catwise-density-levelset-fit.json'))['migkas_fits']
E=json.load(open(N+'catwise-eta-migkas-fit.json'))['migkas_fits']['d14.3']; E['kind']='eta'; E['shape']={'eta':E['eta']}
M=[('Exponential (lock, ℓ₀=15.2°)',F['exp_fixed15.2']),('Soft core (lock, ℓ₀=15.2°)',F['soft_fixed15.2']),('Stretched exp η=0.05 (lock, ℓ₀ free)',E)]
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
fig,axs=plt.subplots(1,3,figsize=(32,13)); rows=[]
for ax,(name,p) in zip(axs,M):
    D=p['D_Gpc']; a=ns['_pc_args'](p['kind'],p['DR'],p['LD'],p.get('shape'))
    rr=pc.run(p['kind'],a,p['mc'],0.0)
    q,rho=bc.geometry(bc.TH1,p['DR'])           # code geometry, units of D
    X,V,VR=pc.integrate(q,rr['mu'],pc.RG,rr['MT'],rr['mub'],rr['Rb'],0.01,200)  # positions today
    Q=q*D; XT=X*D; Rb=rr['Rb']*D; dG=rho*D
    # density map rho/rho(D) from code table S(r), truncated at Rb
    lim=1.9*D if True else 0; g=np.linspace(-lim,lim,700); xx,yy=np.meshgrid(g,g); r=np.hypot(xx,yy)/D
    S=np.interp(r,pc.RG,rr['S']); SD=np.interp(1.0,pc.RG,rr['S']); dens=np.where(r<=rr['Rb'],S/SD,np.nan)
    im=ax.pcolormesh(g,g,dens,cmap='viridis',vmin=np.nanmin(dens),vmax=np.nanmax(dens),shading='auto')
    cb=fig.colorbar(im,ax=ax,fraction=.046); cb.set_label('ρ / ρ(D)  (code profile S(r))',fontsize=18); cb.ax.tick_params(labelsize=15)
    t=np.linspace(0,2*np.pi,400)
    ax.plot(Rb*np.cos(t),Rb*np.sin(t),'w--',lw=2,label=f'profile truncation Rb = {Rb:.1f}')
    ax.plot(0,0,'*',color='r',ms=26,label='CoE (mass profile centre)')
    ax.plot(Q[0,0],Q[0,1],'o',color='orange',ms=16,mec='k',label=f'Observer, D = {D:.2f}')
    ax.plot(Q[1:,0],Q[1:,1],'-',color='w',lw=3,label=f'visible circle (code q), d = {dG:.2f}')
    ax.plot(XT[1:,0],XT[1:,1],':',color='pink',lw=3,label='same points integrated to today')
    ax.plot(XT[0,0],XT[0,1],'x',color='pink',ms=16,mew=3)
    pm=Q[1+180],Q[1+0]   # θ=180 (D−d), θ=0 (D+d)
    ax.plot([pm[0][0],pm[1][0]],[pm[0][1],pm[1][1]],'s',color='magenta',ms=13,mec='k',label=f'D−d = {pm[0][0]:.1f}, D+d = {pm[1][0]:.1f}')
    ax.annotate('',xy=(pm[1][0]+0.35*D,0),xytext=(pm[1][0]+0.05*D,0),arrowprops=dict(arrowstyle='-|>',lw=4,color='yellow'))
    ax.text(pm[1][0]+0.05*D,-0.28*D,f"θ=0: ℓ₀={p['l0']:.1f}°\n(anti-CoE, outward)",color='yellow',fontsize=15)
    ax.annotate('',xy=(Q[0,0]*1.0,0.75*D),xytext=(Q[0,0]*0.6,0.45*D),arrowprops=dict(arrowstyle='-|>',lw=4,color='cyan'))
    ax.text(-0.9*D,0.85*D,'outward from CoE\n(radial launch v ∝ r)',color='cyan',fontsize=15)
    L=p['L_Gpc_locked']
    ax.set_title(f"{name}\nD={D:.2f} Gpc, mc=M(<D)/M_crit={p['mc']:.3f}, L={L:.4g} Gpc",fontsize=19)
    ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim);ax.set_aspect('equal')
    ax.set_xlabel('x (Gpc)',fontsize=19);ax.set_ylabel('y (Gpc)',fontsize=19);ax.tick_params(labelsize=15)
    ax.legend(fontsize=14,loc='lower left',framealpha=.85)
    # verification
    chk=[('CoE',(0,0),(0,0)),('Observer x',D*1.0,Q[0,0]),('D-d x',D*(1-rho),pm[0][0]),('D+d x',D*(1+rho),pm[1][0]),
         ('d (Gpc)',D/p['DR']*ns['RV']/ns['RV']* (ns['RV']/(D/p['DR'])) ,dG),('Rb (Gpc)',1.5*(1+1/p['DR'])*D,Rb),
         ('visible radius max dev',0.0,float(np.max(np.abs(np.hypot(Q[1:,0]-D,Q[1:,1])-dG))))]
    for el,c,v in chk:
        assert np.allclose(c,v,atol=1e-9),(el,c,v); rows.append((p['kind'],el,c,v))
    assert np.isclose(dG,ns['RV'],rtol=1e-6)
fig.suptitle('Fitted universes as the paper code defines them (to scale, Gpc; CoE at origin, θ=0 = ℓ₀)',fontsize=24)
fig.tight_layout();fig.savefig('/workspace/cosmology-course/notes/paper_repro/fitted_universe.png',dpi=110)
print(f"{'model':6s} {'element':24s} {'code':>22s} {'plotted':>22s}")
for k,e,c,v in rows: print(f"{k:6s} {e:24s} {str(np.round(c,6)):>22s} {str(np.round(v,6)):>22s}")
for name,p in M: print(name,'RV',ns['RV'],'Ht',p['Ht'])
