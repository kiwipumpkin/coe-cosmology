import glob,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
O='/workspace/redweb2/realweb/los_flow_R/';ST='/workspace/redweb2/realweb/los/';sm=np.load('/workspace/redweb2/sm_red.npy')
E=sorted(set(np.round(np.concatenate([np.array([1,2,5])*10.**p for p in range(-8,4)]),15)),reverse=True)
l=np.arange(2,31);pm=np.median(Pl[2:,1]);t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
th,ph=hp.pix2ang(32,np.arange(12288));et=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1);ep=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
def grad(m):_,dt,dp=hp.alm2map_der1(hp.map2alm(m,lmax=95),32,lmax=95);return dt,dp
VAR={'grow_full':('grow','mlf'),'fixed_full':('fixed','mlf'),'grow_rel':('grow','mlr'),'fixed_rel':('fixed','mlr')}
GS=grad(sm);_gc={}
def sky(z,R,e,var):
    v,ml=VAR[var];p=f'R{R:g}_';key=(id(z),p+'sw_'+v)
    if key not in _gc:_gc[key]=(grad(z[p+'sw_'+v]),z)
    gw=_gc[key][0];dt=GS[0]+e*gw[0];dp=GS[1]+e*gw[1]
    Tp=sm+e*z[p+'sw_'+v];out=Tp+z[p+'rs_cc']+e*(z[p+'isw_'+v]+z[p+ml+'_'+v])
    al=e*z[p+'al_'+v]+z[p+'ac'];return out+(al*et).sum(1)*dt+(al*ep).sum(1)*dp
def ladder(Ds):
    st=[(np.median(D[4:31]),abs(np.polyfit(l[2:],np.log(np.maximum(D[4:31],1e-30)),1)[0])<0.05) for D in Ds]
    for i in range(len(E)-1):
        if st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0]:
            m1,m2=np.log(st[i][0]),np.log(st[i+1][0]);f=(m1-np.log(pm))/(m1-m2) if m1!=m2 else 0
            return True,float(np.exp(np.log(E[i])+f*(np.log(E[i+1])-np.log(E[i])))),st
    return False,None,st
def lvl(st):
    md=np.array([s[0] for s in st]);i=np.where((md[:-1]>=pm)&(md[1:]<pm))[0]
    if not len(i):return None
    i=i[0];f=(np.log(md[i])-np.log(pm))/(np.log(md[i])-np.log(md[i+1]));return float(np.exp(np.log(E[i])+f*(np.log(E[i+1])-np.log(E[i]))))
pw=lambda m:float(np.sqrt(np.sum(Dl(m)[2:31])))
k='113';z=dict(np.load(O+f'comp/u{k}.npz'));var='grow_full';D0=Dl(sm)
pm5=np.median(Pl[(Pl[:,0]>=5)&(Pl[:,0]<=30),1]);print('Planck median l5-30',pm5)
fig2,ax2=plt.subplots(2,1,figsize=(11,9))
for b,R in zip(ax2,(8,14)):
  Es=np.geomspace(0.01,2,30);rows=[]
  for e in Es:
    D=Dl(sky(z,R,e,var));rows.append((e,np.median(D[5:31]),chi(D),r2(D),D[2:5].copy()))
  for r in rows:print(R,f'E={r[0]:.4g} med={r[1]:.4g} chi2={r[2]:.4g} r2={r[3]:.3f} l2-4={np.round(r[4],1)}')
  md=np.array([r[1] for r in rows]);i=np.where((md[:-1]<pm5)&(md[1:]>=pm5))[0][0];sel=[rows[j] for j in (i+2,i+1,i,i-1)]
  b.plot(l,Lt,'k-',lw=2,label=f'ΛCDM χ²={chi(np.r_[0,0,Lt]):.1f}, r²={r2(np.r_[0,0,Lt]):.2f}');b.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}, r²={r2(D0):.2f}')
  for j,r in enumerate(sel):D=Dl(sky(z,R,r[0],var));b.plot(l,D[2:31],color=plt.cm.viridis(j/3),lw=2,label=f'E={r[0]:.3g} χ²={r[2]:.1f}, r²={r[3]:.2f}');print('SEL',R,r[0],r[1],r[2],r[3],np.round(r[4],1))
  b.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck');b.set_ylabel('D_ℓ (μK²)')
  b.set_ylim(0,2500);b.set_title(f'#{k} {var} R={R} Gpc: 4 closest level sets (linear)');b.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
ax2[-1].set_xlabel('ℓ');fig2.tight_layout();fig2.savefig(O+'flowR_grow_full_lin4_hi.png',dpi=90)
