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
a1=Dl(sky(z,8,0.015,var));_gc.clear();z2=dict(np.load(O+f'comp/u{k}.npz'));a2=Dl(sky(z2,8,0.015,var));print('determinism max|diff|',np.abs(a1-a2).max(),np.array_equal(a1,a2))
pm5=np.median(Pl[(Pl[:,0]>=5)&(Pl[:,0]<=30),1]);Es=np.geomspace(0.02,0.01,8)
fig2,ax2=plt.subplots(3,1,figsize=(11,13.5))
for b,R in zip(ax2,(4,8,14)):
  b.plot(l,Lt,'k-',lw=2,label=f'ΛCDM χ²={chi(np.r_[0,0,Lt]):.1f}, r²={r2(np.r_[0,0,Lt]):.2f}');b.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}, r²={r2(D0):.2f}')
  md=[]
  for j,e in enumerate(Es):
    D=Dl(sky(z,R,e,var));m=np.median(D[5:31]);md.append(m);b.plot(l,D[2:31],color=plt.cm.viridis(j/7),lw=2,label=f'E={e:.3g} χ²={chi(D):.1f}, r²={r2(D):.2f}')
    print(f'R={R} E={e:.4f} med={m:.0f} chi2={chi(D):.1f} r2={r2(D):.3f}')
  md=np.array(md);i=np.where((md[:-1]>=pm5)&(md[1:]<pm5))[0]
  if len(i):i=i[0];f=np.log(md[i]/pm5)/np.log(md[i]/md[i+1]);print('CROSS',R,np.exp(np.log(Es[i])+f*np.log(Es[i+1]/Es[i])))
  else:print('CROSS',R,'none in range')
  b.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck');b.set_ylabel('D_ℓ (μK²)')
  b.set_ylim(0,2500);b.set_title(f'#{k} {var} R={R} Gpc: fine E between 0.01 and 0.02 (linear)');b.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
ax2[-1].set_xlabel('ℓ');fig2.tight_layout();fig2.savefig(O+'flowR_grow_full_lin_fineE.png',dpi=90)
