import glob,json,sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
O='/workspace/redweb2/realweb/los_flow/';ST='/workspace/redweb2/realweb/los/';sm=np.load('/workspace/redweb2/sm_red.npy')
E=list(np.array([1e4,5e3,2e3,1e3,500,200,100,50,20,10,5,2,1,.5,.2,.1,.05,.02,.01,.005,.002,.001,5e-4,2e-4,1e-4,5e-5,2e-5,1e-5,5e-6,2e-6,1e-6,5e-7,2e-7,1e-7]))
l=np.arange(2,31);pm=np.median(Pl[2:,1]);t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
th,ph=hp.pix2ang(32,np.arange(12288));et=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1);ep=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
def grad(m):_,dt,dp=hp.alm2map_der1(hp.map2alm(m,lmax=95),32,lmax=95);return dt,dp
def sky(z,s,e,v,parts=('sw','isw','lens')):
    Tp=sm+(e*z['sw_'+v] if 'sw' in parts else 0);out=Tp+s['rs_cc']
    if 'isw' in parts:out=out+e*(z['isw_'+v]+z['iswdef_'+v])
    if 'lens' in parts:al=e*z['al_'+v]+s['ac'];dt,dp=grad(Tp);out=out+(al*et).sum(1)*dt+(al*ep).sum(1)*dp
    return out
def ladder(Ds):
    st=[(np.median(D[4:31]),abs(np.polyfit(l[2:],np.log(np.maximum(D[4:31],1e-30)),1)[0])<0.05) for D in Ds]
    for i in range(len(E)-1):
        if st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0]:
            m1,m2=np.log(st[i][0]),np.log(st[i+1][0]);f=(m1-np.log(pm))/(m1-m2) if m1!=m2 else 0
            return True,float(np.exp(np.log(E[i])+f*(np.log(E[i+1])-np.log(E[i]))))
    return False,None
SR=json.load(open(ST+'los_results.json'));res={}
VAR=[('grow',('sw','isw','lens'),'grow'),('grow_sw',('sw',),'grow'),('grow_sw_isw',('sw','isw'),'grow'),('grow_sw_lens',('sw','lens'),'grow'),('fixed',('sw','isw','lens'),'fixed')]
for f in sorted(glob.glob(O+'comp/u*.npz')):
    k=f.split('/u')[1][:-4];z=dict(np.load(f));s=dict(np.load(ST+f'comp/u{k}.npz'));r={}
    for nm,p,v in VAR:
        D=[Dl(sky(z,s,e,v,p)) for e in E];ok,es=ladder(D);r[nm]=dict(guess=ok,E=es)
        if nm in('grow','fixed'):
            c2=[chi(x) for x in D];jb=int(np.argmin(c2));r[nm].update(chi2_best=c2[jb],r2_best=r2(D[jb]),E_best=float(E[jb]))
            if es:Dx=Dl(sky(z,s,es,v,p));r[nm].update(chi2_Estar=chi(Dx),r2_Estar=r2(Dx))
    # rms of each component (dipole-removed, l>=2 power) at E=1
    for v in ('grow','fixed'):r['rms_'+v]={q:float(np.sqrt(np.sum(Dl(z[q+'_'+v])[2:31]))) for q in ('sw','isw','iswdef')}
    r['static']=SR[k]['los'];r.update({q:SR[k][q] for q in ('ss4','ss8','ss14')});res[k]=r;print(k,{n:r[n]['E'] for n,_,_ in VAR},flush=True)
json.dump(res,open(O+'flow_results.json','w'))
keys=['grow','grow_sw','grow_sw_isw','grow_sw_lens','fixed','static','ss4','ss8','ss14'];summ={}
def st_(v):return dict(n_guess=len(v),median_E=float(np.median(v)) if v else None,p16=float(np.percentile(v,16)) if v else None,p84=float(np.percentile(v,84)) if v else None)
for q in keys:summ[q]=st_([r[q]['E'] for r in res.values() if r[q]['guess']]);summ[q]['N']=len(res)
for v in ('grow','fixed'):
    for q in ('chi2_Estar','r2_Estar','chi2_best','r2_best'):summ[v][q+'_median']=float(np.median([r[v][q] for r in res.values() if q in r[v]]))
    summ['rms_'+v]={q:float(np.median([r['rms_'+v][q] for r in res.values()])) for q in ('sw','isw','iswdef')}
json.dump(summ,open(O+'flow_summary.json','w'),indent=1);print(json.dumps(summ,indent=1))
lab=['flow grow\n(all)','grow\nSW','grow\nSW+ISW','grow\nSW+lens','flow fixed-\nmass (all)','static\nlight path','single\nsphere 4','single\nsphere 8','single\nsphere 14']
fig,ax=plt.subplots(1,2,figsize=(16,6))
for i,q in enumerate(keys):
    v=[r[q]['E'] for r in res.values() if r[q]['guess']]
    ax[0].scatter(np.full(len(v),i)+np.random.uniform(-.15,.15,len(v)),v,s=10,alpha=.6)
    if v:ax[0].errorbar(i,summ[q]['median_E'],yerr=[[summ[q]['median_E']-summ[q]['p16']],[summ[q]['p84']-summ[q]['median_E']]],fmt='k_',ms=25,capsize=6);ax[0].text(i+.2,summ[q]['median_E'],f"{summ[q]['median_E']:.3g}",fontsize=8)
ax[0].set_yscale('log');ax[0].set_xticks(range(len(keys)),lab,fontsize=8);ax[0].set_ylabel('required strength E* (bracket interp.)');ax[0].set_title(f'Required web strength, {len(res)} MC universes (median, 16–84%)')
ax[1].bar(range(len(keys)),[summ[q]['n_guess'] for q in keys]);ax[1].set_xticks(range(len(keys)),lab,fontsize=8);ax[1].set_ylabel(f'bracket-rule guess count / {len(res)}')
for i,q in enumerate(keys):ax[1].text(i,summ[q]['n_guess'],str(summ[q]['n_guess']),ha='center',va='bottom')
ax[1].set_title('Bracket-rule guesses (count)');plt.tight_layout();plt.savefig(O+'flow_vs_static_strength.png',dpi=100)
D0=Dl(sm)
for v in ('grow','fixed'):
  ks=list(res)[:4];fig,ax=plt.subplots(4,1,figsize=(11,20));fig2,ax2=plt.subplots(4,1,figsize=(11,18))
  for a,b,k in zip(ax,ax2,ks):
    z=dict(np.load(O+f'comp/u{k}.npz'));s=dict(np.load(ST+f'comp/u{k}.npz'));es=res[k][v]['E'] or res[k][v]['E_best']
    jc=int(np.argmin(abs(np.log(np.array(E))-np.log(es))));sub=E[max(0,jc-5):jc+6]
    for A in (a,b):A.plot(l,Lt,'k-',lw=2,label=f'ΛCDM χ²={chi(np.r_[0,0,Lt]):.1f}, r²={r2(np.r_[0,0,Lt]):.2f}');A.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}, r²={r2(D0):.2f}')
    for j,e in enumerate(sub):D=Dl(sky(z,s,e,v));a.plot(l,D[2:31],color=plt.cm.viridis(j/(len(sub)-1)),label=f'E={e:.3g} χ²={chi(D):.4g}, r²={r2(D):.2f}')
    fe=es*np.geomspace(.5,2,13);FD=sorted([(chi(Dl(sky(z,s,e,v))),e) for e in fe])[:4]
    for j,(c2,e) in enumerate(sorted(FD,key=lambda x:-x[1])):D=Dl(sky(z,s,e,v));b.plot(l,D[2:31],color=plt.cm.viridis(j/3),lw=2,label=f'E={e:.3g} χ²={c2:.1f}, r²={r2(D):.2f}')
    for A in (a,b):A.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck');A.set_ylabel('D_ℓ (μK²)')
    g=res[k][v];a.set_yscale('log');a.set_ylim(1e-3,None);a.set_title(f'web riding CoE flow ({v}) #{k}: bracket guess {"in bracket, E*≈%.3g"%g["E"] if g["guess"] else "no bracket"}')
    a.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85);b.set_ylim(0,2500);b.set_title(f'#{k} flow ({v}), 4 closest level sets (linear, Planck range)');b.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
  ax[-1].set_xlabel('ℓ');ax2[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig(O+f'flow_{v}_ladder.png',dpi=90);fig2.tight_layout();fig2.savefig(O+f'flow_{v}_lin4.png',dpi=90);plt.close('all')
