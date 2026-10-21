import glob,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
O='/workspace/redweb2/realweb/los/';sm=np.load('/workspace/redweb2/sm_red.npy');S=np.load('/workspace/redweb2/realweb/mc/mc_sky.npy',mmap_mode='r')
E=[100,50,20,10,5,2,1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];l=np.arange(2,31);pm=np.median(Pl[2:,1])
t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
th,ph=hp.pix2ang(32,np.arange(12288));et=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1);ep=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
def grad(m):_,dt,dp=hp.alm2map_der1(hp.map2alm(m,lmax=95),32,lmax=95);return dt,dp
def sky(z,e,parts=('sw','rs','lens')):
    Tp=sm+(e*z['sw'] if 'sw' in parts else 0);out=Tp.copy()
    if 'rs' in parts:out=out+e*e*z['rs_ww']+e*z['rs_wc']+z['rs_cc']
    if 'lens' in parts:
        al=e*z['aw']+z['ac'];dt,dp=grad(Tp);out=out+(al*et).sum(1)*dt+(al*ep).sum(1)*dp
    return out
def ladder(Ds):
    st=[(np.median(D[4:31]),abs(np.polyfit(l[2:],np.log(np.maximum(D[4:31],1e-30)),1)[0])<0.05) for D in Ds]
    for i in range(len(E)-1):
        if st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0]:
            m1,m2=np.log(st[i][0]),np.log(st[i+1][0]);f=(m1-np.log(pm))/(m1-m2) if m1!=m2 else 0
            return True,float(np.exp(np.log(E[i])+f*(np.log(E[i+1])-np.log(E[i]))))
    return False,None
res={};files=sorted(glob.glob(O+'comp/u*.npz'))
for f in files:
    k=int(f.split('/u')[1][:-4]);z=dict(np.load(f));D=[Dl(sky(z,e)) for e in E];ok,es=ladder(D)
    r=dict(los=dict(guess=ok,E=es,chi2=[chi(x) for x in D],r2=[r2(x) for x in D]),floor=chi(Dl(sky(z,0))))
    for p,nm in ((('sw',),'sw_only'),(('sw','rs'),'sw_rs'),(('sw','lens'),'sw_lens')):
        o2,e2=ladder([Dl(sky(z,e,p)) for e in E]);r[nm]=dict(guess=o2,E=e2)
    EE=E;E=E[6:]
    for iR,R in ((3,4.0),(4,8.0),(5,14.0)):
        o2,e2=ladder([Dl(sm+e*np.array(S[k][iR])) for e in E]);r[f'ss{R:g}']=dict(guess=o2,E=e2)
    E=EE
    res[k]=r;print(k,ok,es,r['ss4']['E'],r['ss8']['E'],r['ss14']['E'],flush=True)
json.dump(res,open(O+'los_results.json','w'))
# summary
keys=['los','sw_only','sw_rs','sw_lens','ss4','ss8','ss14'];summ={}
for q in keys:
    v=[r[q]['E'] for r in res.values() if r[q]['guess']];summ[q]=dict(n_guess=len(v),N=len(res),median_E=float(np.median(v)) if v else None,p16=float(np.percentile(v,16)) if v else None,p84=float(np.percentile(v,84)) if v else None)
summ['floor_chi2_median']=float(np.median([r['floor'] for r in res.values()]))
json.dump(summ,open(O+'los_summary.json','w'),indent=1);print(json.dumps(summ,indent=1))
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
lab=['light path\n(all)','SW only','SW+RS','SW+lens','single\nsphere 4','single\nsphere 8','single\nsphere 14']
for i,q in enumerate(keys):
    v=[r[q]['E'] for r in res.values() if r[q]['guess']]
    ax[0].scatter(np.full(len(v),i)+np.random.uniform(-.15,.15,len(v)),v,s=10,alpha=.6)
    if v:ax[0].errorbar(i,summ[q]['median_E'],yerr=[[summ[q]['median_E']-summ[q]['p16']],[summ[q]['p84']-summ[q]['median_E']]],fmt='k_',ms=25,capsize=6)
ax[0].set_yscale('log');ax[0].set_xticks(range(len(keys)),lab,fontsize=8);ax[0].set_ylabel('required strength E* (bracket interp.)');ax[0].set_title(f'Required web strength, {len(res)} MC universes (median, 16–84%)')
ax[1].bar(range(len(keys)),[summ[q]['n_guess'] for q in keys]);ax[1].set_xticks(range(len(keys)),lab,fontsize=8);ax[1].set_ylabel(f'bracket-rule guess count / {len(res)}')
for i,q in enumerate(keys):ax[1].text(i,summ[q]['n_guess'],str(summ[q]['n_guess']),ha='center',va='bottom')
ax[1].set_title('Bracket-rule guesses (count)');plt.tight_layout();plt.savefig(O+'los_vs_single_strength.png',dpi=100)
# ladders + linear 4-closest for 4 random universes
ks=list(res)[:4];D0=Dl(sm)
fig,ax=plt.subplots(len(ks),1,figsize=(11,5*len(ks)))
fig2,ax2=plt.subplots(len(ks),1,figsize=(11,4.5*len(ks)))
for a,b,k in zip(ax,ax2,ks):
    z=dict(np.load(O+f'comp/u{k}.npz'));Ds=[Dl(sky(z,e)) for e in E];Df=Dl(sky(z,0))
    for A in (a,b):A.plot(l,Lt,'k-',lw=2,label=f'ΛCDM χ²={chi(np.r_[0,0,Lt]):.1f}, r²={r2(np.r_[0,0,Lt]):.2f}');A.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}, r²={r2(D0):.2f}')
    a.plot(l,Df[2:31],color='orange',ls=':',lw=2,label=f'red + CoE light path, no web χ²={chi(Df):.3g}, r²={r2(Df):.2f}')
    for j,(e,D) in enumerate(zip(E,Ds)):a.plot(l,D[2:31],color=plt.cm.viridis(j/(len(E)-1)),label=f'{e*100:g}% χ²={chi(D):.4g}, r²={r2(D):.2f}')
    # fine ladder around best, keep 4 closest
    jb=int(np.argmin([chi(D) for D in Ds]));fe=E[jb]*np.geomspace(.4,2.1,12);FD=[(chi(Dl(sky(z,e))),e) for e in fe];FD=sorted(FD)[:4]
    for j,(c2,e) in enumerate(sorted(FD,key=lambda x:-x[1])):D=Dl(sky(z,e));b.plot(l,D[2:31],color=plt.cm.viridis(j/3),lw=2,label=f'{e*100:.3g}% χ²={c2:.1f}, r²={r2(D):.2f}')
    for A in (a,b):A.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck');A.set_ylabel('D_ℓ (μK²)')
    g=res[k]['los'];tt=f'light path: red + MC web #{k}  (bracket guess: {"in bracket" if g["guess"] else "no bracket"}' + (f', E*≈{g["E"]*100:.2g}%)' if g['E'] else ')')
    a.set_yscale('log');a.set_ylim(1e-3,None);a.set_title(tt);a.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
    b.set_ylim(0,2500);b.set_title(f'web #{k} light path, 4 closest level sets (linear, Planck range)');b.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
ax[-1].set_xlabel('ℓ');ax2[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig(O+'los_ladder.png',dpi=90);fig2.tight_layout();fig2.savefig(O+'los_lin4.png',dpi=90)
