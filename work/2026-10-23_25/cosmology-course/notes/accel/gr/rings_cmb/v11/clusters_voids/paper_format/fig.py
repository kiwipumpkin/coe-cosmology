import os,sys,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.interpolate import CubicSpline
sys.path.insert(0,'../../rotfit');import fit
SEQ='/workspace/clump_scan/sequential';c=os.getcwd();os.chdir(SEQ);sys.path.insert(0,SEQ);import seq_core as S;os.chdir(c)
M3='../leftover_vs_n/m345/';R=json.load(open(M3+'res_both_m345.json'));SUM=json.load(open(M3+'m345_summary.json'))
TYP={'mass':dict(dl=30,dh=300,rmax=0.30),'void':dict(dl=-0.45,dh=-0.25,rmax=2.14)};DL=14.3
nside=512;V=np.array(hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))).T
o=np.load('../../shrink_rec/b_1100_k0.npz')['o'];aC=fit.vec(0,15);CoE=-aC
m=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(np.degrees(np.arccos(np.clip(V@aC,-1,1)))));m=m/m.mean()-1
B=[]
for p,t in zip(R['params'],R['types']):
    T=TYP[t];al=p[2]*T['rmax']/DL;d=T['dl']+p[3]*(T['dh']-T['dl']);n=hp.ang2vec(p[0],p[1]);B.append(dict(t=t,n=n,al=al,amp=np.log1p(d)/3,d=d,r=p[2]*T['rmax']))
ma=m.copy()
for b in B:
    pix=hp.query_disc(nside,b['n'],b['al'],inclusive=False);ma[pix]+=b['amp']
def strip(x):return hp.remove_dipole(x,fitval=False,verbose=False) if 'verbose' in hp.remove_dipole.__code__.co_varnames else hp.remove_dipole(x)
mb,mA=strip(m),strip(ma)
def lr(x):cl=hp.anafast(x,lmax=12);l=np.arange(13);return np.sqrt((2*l+1)*cl/(4*np.pi))/1e-5
Lb,La=lr(mb),lr(mA);g=lambda L,ls:float(np.sqrt(np.sum(L[ls]**2)))
num=dict(before=dict(l2=float(Lb[2]),l345=g(Lb,[3,4,5]),l6_12=g(Lb,list(range(6,13))),total_sky_rms=float(np.std(mb)/1e-5)),after=dict(l2=float(La[2]),l345=g(La,[3,4,5]),l6_12=g(La,list(range(6,13))),total_sky_rms=float(np.std(mA)/1e-5)))
# ring: great circle through CoE, perpendicular toward galactic north
w=np.array([0,0,1.])-(np.array([0,0,1.])@CoE)*CoE;w/=np.linalg.norm(w);ph=np.radians(np.arange(0,360,0.1))
Nr=np.outer(np.cos(ph),CoE)+np.outer(np.sin(ph),w);pr=hp.vec2pix(nside,*Nr.T);rb,ra=mb[pr]/1e-5,mA[pr]/1e-5
num['ring']=dict(before_rms=float(np.std(rb)),after_rms=float(np.std(ra)),before_range=[float(rb.min()),float(rb.max())],after_range=[float(ra.min()),float(ra.max())])
hits=[i for i,b in enumerate(B) if np.min(np.degrees(np.arccos(np.clip(Nr@b['n'],-1,1))))<np.degrees(b['al'])];num['bodies_on_ring_circle']=hits
json.dump(num,open('before_after_numbers.json','w'),indent=1);print(num)
# web (seed 3), Chris's frame rotated so CoE is left
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X,r,M=f['Xp'],f['rp'],f['Mp'];isV=f['kp']=='void';CG=np.array(S.C_G);A=float(S.A_IN_G)
u=CG/np.linalg.norm(CG);wz=np.cross(u,[0,0,1.]);wz/=np.linalg.norm(wz);zz=np.cross(u,wz);P=np.c_[-(X@u),X@wz,X@zz];sl=np.abs(P[:,2])<3
# our added bodies in our frame (x along CoE mirrored, y along w)
zb=np.cross(CoE,w);BP=np.array([[-(DL*b['n'])@CoE,(DL*b['n'])@w,(DL*b['n'])@zb] for b in B])
def slab(ax,after):
    ax.scatter(P[sl&isV,0],P[sl&isV,1],s=10,c='#2166ac',label='void centre');ax.scatter(P[sl&~isV,0],P[sl&~isV,1],s=5*M[sl&~isV]/np.median(M[~isV]),c='#d62728',label='cluster')
    ax.add_patch(Circle((0,0),14.3,fill=False,ec='k',lw=2.2));ax.add_patch(Circle((-np.linalg.norm(CG),0),A,fc='pink',ec='m',alpha=.55))
    ax.annotate('',xy=(-30,0),xytext=(0,0),arrowprops=dict(arrowstyle='-|>',color='#a93226',lw=2.5,mutation_scale=22));ax.text(-30,1.3,'to CoE',color='#a93226',fontsize=14,fontweight='bold');ax.plot(0,0,'k*',ms=14)
    if after:
        for i,(b,q) in enumerate(zip(B,BP)):
            ins=abs(q[2])<3;col='#d62728' if b['t']=='mass' else '#2166ac'
            ax.add_patch(Circle(q[:2],max(b['r'],0.45),fc=col if ins else 'none',ec='k',lw=2.5,zorder=6));ax.text(q[0]+0.8,q[1]+0.8,f"{i+1}",fontsize=13,fontweight='bold',zorder=7)
    ax.set_xlim(-32,32);ax.set_ylim(-32,32);ax.set_aspect('equal');ax.tick_params(labelsize=12);ax.set_xlabel('Gpc (CoE to the left)',fontsize=13)
fig=plt.figure(figsize=(20,21));g1=fig.add_gridspec(2,2,height_ratios=[1.25,0.85],hspace=0.22)
a1=fig.add_subplot(g1[0,0]);slab(a1,False);a1.set_title('Before: random compensated web, seed 3 (no body touches the ring)',fontsize=15)
a2=fig.add_subplot(g1[0,1]);slab(a2,True);a2.set_title('After: + 8 bodies placed on the ring (black outline, numbered;\nfilled = within the 3 Gpc slab, open = projected from off the slab)',fontsize=15)
h=[plt.Line2D([],[],marker='o',ls='',c='#2166ac'),plt.Line2D([],[],marker='o',ls='',c='#d62728'),plt.Line2D([],[],c='k',lw=2),plt.Line2D([],[],marker='s',ls='',mfc='pink',mec='m',ms=14),plt.Line2D([],[],marker='o',ls='',mfc='w',mec='k',mew=2.5,ms=14)]
fig.legend(h,['void centre','cluster (size ∝ mass)','our sphere / CMB ring, 14.3 Gpc','CoE core cap','added body (6 masses, 2 voids)'],loc='upper center',ncol=5,fontsize=14,frameon=False,bbox_to_anchor=(.5,.93))
a3=fig.add_subplot(g1[1,:]);x=np.degrees(ph)
a3.axhspan(-1,1,color='#b2f2d4',alpha=.8,label='±1× (10⁻⁵)');a3.plot(x,rb,color='.55',lw=3,label=f"before: ring RMS {num['ring']['before_rms']:.0f}×")
a3.plot(x,ra,color='#08306b',lw=1.8,label=f"after: ring RMS {num['ring']['after_rms']:.0f}× (bodies cut this circle: {len(hits)} of 8)")
a3.set_xlim(0,360);a3.set_xticks(range(0,361,30));a3.set_xlabel('angle around the ring from the CoE direction (°)',fontsize=15);a3.set_ylabel('residual (units of 10⁻⁵), dipole removed by motion',fontsize=14);a3.tick_params(labelsize=13);a3.grid(alpha=.3)
a3.legend(loc='lower center',bbox_to_anchor=(.5,1.0),ncol=3,fontsize=14,frameon=False)
b,a=num['before'],num['after']
cap=(f"Our model (mass b, CatWISE slope, 1+z ≈ 1100), response dT/T = ln(1+δ)/3 where a body cuts the recombination surface at 14.3 Gpc. Dipole absorbed by our motion (637 km/s) in both curves.\n"
 f"Full sky, units of 10⁻⁵: before l=2 {b['l2']:.0f}, l=3–5 {b['l345']:.1f}, l=6–12 {b['l6_12']:.2f};  after l=2 {a['l2']:.0f} (not targeted, still dominant), l=3–5 {a['l345']:.2f} (cancelled), l=6–12 {a['l6_12']:.0f} (raised by the bodies).\n"
 "The seed-3 web contributes nothing under our rule: its generator keeps every body off the ring. Ring = great circle through the CoE; slab frames of the web and of the added bodies are each aligned with CoE to the left (perpendicular orientation not matched).")
fig.text(.5,.06,cap,ha='center',va='top',fontsize=12.5)
fig.savefig('ring_residual_before_after.png',dpi=70,bbox_inches='tight');plt.close(fig)
# FIG 71 style: leftover alone + mode amplitudes
fig,axs=plt.subplots(1,2,figsize=(22,7.5),gridspec_kw=dict(width_ratios=[2,1]))
axs[0].axhspan(-1,1,color='#b2f2d4',alpha=.8,label='±1×');axs[0].plot(x,ra,color='#08306b',lw=1.8,label='leftover after the 8 bodies');axs[0].set_xlim(0,360);axs[0].set_xticks(range(0,361,30))
axs[0].set_xlabel('angle around the ring from the CoE (°)',fontsize=14);axs[0].set_ylabel('residual (10⁻⁵)',fontsize=14);axs[0].legend(loc='lower center',bbox_to_anchor=(.5,1.0),ncol=2,fontsize=13,frameon=False);axs[0].grid(alpha=.3)
l=np.arange(2,13);axs[1].bar(l-0.2,Lb[2:],0.4,color='.6',label='before');axs[1].bar(l+0.2,La[2:],0.4,color='#08306b',label='after');axs[1].set_yscale('log');axs[1].axhline(1,color='g',ls='--',label='1×')
axs[1].set_xticks(l);axs[1].set_xlabel('multipole l (full sky)',fontsize=14);axs[1].set_ylabel('RMS per l (10⁻⁵)',fontsize=14);axs[1].legend(loc='lower center',bbox_to_anchor=(.5,1.0),ncol=3,fontsize=13,frameon=False)
fig.text(.5,-.06,f"Leftover shown alone (as in Chris's Fig. 71). l = 3–5 cancelled to {a['l345']:.2f}×; the quadrupole ({a['l2']:.0f}×) still dominates, and l = 6–12 ({a['l6_12']:.0f}×) is raised by the added bodies.",ha='center',fontsize=13)
fig.savefig('leftover_alone.png',dpi=80,bbox_inches='tight')
