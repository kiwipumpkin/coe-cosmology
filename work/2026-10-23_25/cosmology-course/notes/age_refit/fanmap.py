import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
law,mode,vm=sys.argv[1:4]; os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
sys.argv=['x',law,mode]; exec(open('fit_lc.py').read().split("ck=f'fwd")[0]); import grltb as G
x=np.array(json.load(open(f'gr/ck_{law}_{mode}_{vm}.json'))['x']); P,o=build(x); f=G.fields(o); L_=F.local(o); lam=L_['Hbar']/H70
io=o['io']; D=o['r'][io]; trec=np.interp(D/1100.,o['rs'][:,io],o['ts'])*F.C
tr=G.trace(f,np.arange(0,360,2.),nstep=500); dl=f['tn']-trec
back_ok=[bool(dl<=tr['lb'][-1] and tr['alive'][min(np.searchsorted(tr['lb'],dl),len(tr['lb'])-1),i]) for i in range(180)]
fc=G.forward_check(f,tr,trec,nstep=1200) if dl<=tr['lb'][-1] else None
fm=G.fan_forward(f,trec,nlab=200,npsi=91,nstep=800)
H=fm['hits'] if fm else []
out=dict(tag=f'{law}_{mode}_{vm}',trec_Gyr=trec/F.C*lam,n_backward_rays_reaching=int(sum(back_ok)),n_forward_hits=len(H),
  theta_covered=sorted(set(int(h['theta'])//10*10 for h in H)),hits=H,fwd_check=({k:np.asarray(v).tolist() for k,v in fc.items()} if fc else None))
json.dump(out,open(f'gr/fan_{law}_{mode}_{vm}.json','w'),indent=1,default=float)
fig,ax=plt.subplots(2,2,figsize=(16,10)); ax=ax.ravel(); thb=np.arange(0,360,2.)
if fc:
    ax[0].semilogy(thb,fc['emit_label'],'o',ms=5,label='backward: emitter label at E_rec'); ax[1].plot(thb,fc['emit_phi'],'o',ms=5,label='backward: emission φ')
    ax[2].semilogy(thb,fc['z_back'],'o',ms=5,label='z backward'); ax[2].semilogy(thb,np.abs(fc['z_fwd']),'x',ms=7,label='|z| forward re-run')
    ax[3].semilogy(thb,np.abs(fc['label_miss'])+1e-12,'o',ms=5,label='forward re-run: |label miss| at today')
for a in ax: a.legend(fontsize=12)
ax=ax[:3]
if H:
    th=np.array([h['theta'] for h in H]); ax[0].plot(th,[h['lab'] for h in H],'o',ms=6); ax[1].plot(th,[h['emit_phi'] for h in H],'o',ms=6); ax[2].semilogy(th,[max(h['z'],1e-9) for h in H],'o',ms=6)
th_b=np.arange(0,360,2.)[np.array(back_ok)]; 
for a,yl in zip(ax,['emitter label r0','emission φ (rad)','z']): a.set_xlabel('arrival sky angle θ (0=away CoE)',fontsize=16);a.set_ylabel(yl,fontsize=16);a.set_xlim(0,360);a.tick_params(labelsize=14);a.grid()
fig.suptitle(f"Forward fan from E_rec ({out['tag']}): {len(H)} emitters reach us today; backward rays reaching E_rec: {sum(back_ok)}/180",fontsize=16)
fig.tight_layout();fig.savefig(f'gr/fan_{law}_{mode}_{vm}.png',dpi=80)
print(out['tag'],'trec',out['trec_Gyr'],'hits',len(H),'back',sum(back_ok),'theta bins',out['theta_covered'])
