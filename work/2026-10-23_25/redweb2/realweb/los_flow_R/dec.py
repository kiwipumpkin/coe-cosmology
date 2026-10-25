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
z=dict(np.load(O+'comp/u113.npz'));var='grow_full';LM=8;ls=np.arange(LM+1)
def comps(R,e):
    p=f'R{R:g}_';gw=grad(z[p+'sw_grow']);dt=GS[0]+e*gw[0];dp=GS[1]+e*gw[1];al=e*z[p+'al_grow']+z[p+'ac']
    C={'smooth red':sm,'SW (web)':e*z[p+'sw_grow'],'RS/cc':z[p+'rs_cc'],'ISW':e*z[p+'isw_grow'],'moving lens':e*z[p+'mlf_grow'],'lensing/bend':(al*et).sum(1)*dt+(al*ep).sum(1)*dp}
    return C
fac=ls*(ls+1)*T0**2/(2*np.pi)
th0,ph0=np.radians(90+15),np.radians(180)   # CoE (l180,b-15)
def coe(a):b=a.copy();hp.rotate_alm(b,0,th0,0);return b   # tested below
# test rotation: map = cos(angle to CoE) -> should give only m=0
v0=hp.ang2vec(th0,ph0);tst=hp.map2alm(hp.pix2vec(32,np.arange(12288))[0]*0+np.dot(np.array(hp.pix2vec(32,np.arange(12288))).T,v0),lmax=LM);tr=coe(tst)
print('rot test |a10|,|a11|:',abs(tr[hp.Alm.getidx(LM,1,0)]),abs(tr[hp.Alm.getidx(LM,1,1)]))
def mpow(a,L):return np.array([abs(a[hp.Alm.getidx(LM,L,m)])**2*(1 if m==0 else 2) for m in range(L+1)])
fig,axs=plt.subplots(3,1,figsize=(12,13));out={}
for ax,(R,e) in zip(axs,((4,0.0135),(8,0.025),(14,0.032))):
    C=comps(R,e);A={k:hp.map2alm(np.asarray(hp.remove_dipole(v)),lmax=LM) for k,v in C.items()}
    tot=sum(A.values());Dt=fac*hp.alm2cl(tot);print(f'\n=== R={R} E={e} total D_l l2-8:',np.round(Dt[2:],1),' check Dl:',np.round(Dl(sum(C.values()))[2:9],1))
    ks=list(A);auto={k:fac*hp.alm2cl(A[k]) for k in ks};cross={}
    for i in range(len(ks)):
        for j in range(i+1,len(ks)):cross[ks[i]+' x '+ks[j]]=2*fac*hp.alm2cl(A[ks[i]],A[ks[j]])
    for k in ks:print(f'{k:14s}',np.round(auto[k][2:],1))
    big={k:v for k,v in cross.items() if np.abs(v[2:]).max()>1}
    for k,v in big.items():print(f'X {k:28s}',np.round(v[2:],1))
    print('l5 Planck ~',Pl[Pl[:,0]==5,1])
    for k in ks:
        r=coe(A[k]);mp=mpow(r,5);print(f' CoE-frame m-fraction l5 {k:14s}',np.round(mp/mp.sum(),3) if mp.sum()>0 else 0)
    r=coe(tot);mp=mpow(r,5);print(' CoE-frame m-fraction l5 TOTAL',np.round(mp/mp.sum(),3))
    rs=coe(A['smooth red']);print(' smooth red CoE-frame m=0 fraction per l:',[round(mpow(rs,L)[0]/mpow(rs,L).sum(),3) for L in range(2,9)],' |a_l0|^2 (D-units):',np.round([fac[L]*abs(rs[hp.Alm.getidx(LM,L,0)])**2/(2*L+1) for L in range(2,9)],1))
    items=[(k,auto[k]) for k in ks if np.abs(auto[k][2:]).max()>0.5]+[('X '+k,v) for k,v in big.items()]
    w=0.8/len(items)
    for n,(k,v) in enumerate(items):ax.bar(np.arange(2,9)+(n-len(items)/2)*w,v[2:],w,label=k)
    ax.plot(np.arange(2,9),Dt[2:],'k*',ms=12,label='total');pp=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=8)];ax.plot(pp[:,0],pp[:,1],'ro',label='Planck')
    ax.set_yscale('symlog',linthresh=10);ax.axhline(0,color='k',lw=.5);ax.set_title(f'#113 grow_full R={R} Gpc E={e}: D_ℓ contributions (auto + 2·cross)');ax.set_ylabel('D_ℓ (μK²)');ax.legend(loc='lower center',ncol=4,fontsize=7,framealpha=.85)
axs[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig(O+'l5_decomp.png',dpi=90)
