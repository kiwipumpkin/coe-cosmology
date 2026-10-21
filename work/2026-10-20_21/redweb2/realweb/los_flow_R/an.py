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
SR=json.load(open(ST+'los_results.json'))
def one(f):
    k=f.split('/u')[1][:-4];z=dict(np.load(f));r={}
    for R in (4,8,14):
        p=f'R{R:g}_';r[R]={'rms':{q:pw(z[p+q]) for q in ('sw_grow','isw_grow','mlf_grow','mlr_grow','sw_fixed','isw_fixed','mlf_fixed','mlr_fixed','dop_src','rs_cc')},'vr':[float(z[p+'vr_min']),float(z[p+'vr_max'])]}
        for var in VAR:
            D=[Dl(sky(z,R,e,var)) for e in E];ok,es,st=ladder(D);el=lvl(st);d=dict(guess=ok,E=es,E_level=el)
            eu=es or el
            if eu:Dx=Dl(sky(z,R,eu,var));d.update(chi2=chi(Dx),r2=r2(Dx),slope=float(np.polyfit(l[2:],np.log(Dx[4:31]),1)[0]))
            r[R][var]=d
    r['static']=SR[k]['los'];r.update({q:SR[k][q] for q in ('ss4','ss8','ss14')});print(k,flush=True);return k,r
from multiprocessing import Pool
with Pool(8) as P_:res=dict(P_.map(one,sorted(glob.glob(O+'comp/u*.npz'))))
res={k:{(int(a) if str(a).isdigit() else a):b for a,b in r.items()} for k,r in res.items()}
json.dump(res,open(O+'flowR_results.json','w'))
def st_(v):return dict(n=len(v),median=float(np.median(v)),p16=float(np.percentile(v,16)),p84=float(np.percentile(v,84))) if v else dict(n=0)
summ={}
for R in (4,8,14):
    for var in VAR:
        rr=[r[R][var] for r in res.values()]
        summ[f'R{R}_{var}']=dict(n_guess=sum(x['guess'] for x in rr),Estar=st_([x['E'] for x in rr if x['guess']]),E_level=st_([x['E_level'] for x in rr if x['E_level']]),
          chi2_median=float(np.median([x['chi2'] for x in rr if 'chi2' in x])),r2_median=float(np.median([x['r2'] for x in rr if 'r2' in x])),slope_median=float(np.median([x['slope'] for x in rr if 'slope' in x])))
    summ[f'R{R}_rms_E1']={q:float(np.median([r[R]['rms'][q] for r in res.values()])) for q in res[list(res)[0]][R]['rms']}
    summ[f'R{R}_vr']=res[list(res)[0]][R]['vr']
for q in ('static','ss4','ss8','ss14'):summ[q]=st_([r[q]['E'] for r in res.values() if r[q]['guess']])
json.dump(summ,open(O+'flowR_summary.json','w'),indent=1);print(json.dumps(summ,indent=1))
# strength chart (E* if guessed else level-match E, marked)
fig,ax=plt.subplots(figsize=(17,6));labs=[];i=0
for R in (4,8,14):
  for var in VAR:
    rr=[r[R][var] for r in res.values()];v=[x['E'] for x in rr if x['guess']];vl=[x['E_level'] for x in rr if x['E_level'] and not x['guess']]
    ax.scatter(np.full(len(v),i)+np.random.uniform(-.15,.15,len(v)),v,s=10,c='C0');ax.scatter(np.full(len(vl),i)+np.random.uniform(-.15,.15,len(vl)),vl,s=10,c='C1',marker='x')
    labs.append(f'R{R} {var}\nguess {len(v)}/{len(rr)}');i+=1
for q,nm in (('static','static LOS'),('ss4','single 4'),('ss8','single 8'),('ss14','single 14')):
    v=[r[q]['E'] for r in res.values() if r[q]['guess']];ax.scatter(np.full(len(v),i)+np.random.uniform(-.15,.15,len(v)),v,s=10,c='C0');labs.append(f'{nm}\nguess {len(v)}/{len(res)}');i+=1
ax.set_yscale('log');ax.set_xticks(range(len(labs)),labs,fontsize=7);ax.set_ylabel('strength E');ax.set_title('blue dots: bracket E* (guess);  orange x: level-match E (no bracket, median ℓ4–30 = Planck median)')
plt.tight_layout();plt.savefig(O+'flowR_strength.png',dpi=100);plt.close()
D0=Dl(sm);k=list(res)[0];z=dict(np.load(O+f'comp/u{k}.npz'))
for var in VAR:
  fig,ax=plt.subplots(3,1,figsize=(11,15));fig2,ax2=plt.subplots(3,1,figsize=(11,13.5))
  for a,b,R in zip(ax,ax2,(4,8,14)):
    g=res[k][R][var];ec=g['E'] or g['E_level'];jc=int(np.argmin(abs(np.log(np.array(E))-np.log(ec))));sub=E[max(0,jc-5):jc+6]
    for A in (a,b):A.plot(l,Lt,'k-',lw=2,label=f'ΛCDM χ²={chi(np.r_[0,0,Lt]):.1f}, r²={r2(np.r_[0,0,Lt]):.2f}');A.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}, r²={r2(D0):.2f}')
    for j,e in enumerate(sub):D=Dl(sky(z,R,e,var));a.plot(l,D[2:31],color=plt.cm.viridis(j/(len(sub)-1)),label=f'E={e:.3g} χ²={chi(D):.4g}, r²={r2(D):.2f}')
    FD=sorted([(chi(Dl(sky(z,R,e,var))),e) for e in ec*np.geomspace(.5,2,13)])[:4]
    for j,(c2,e) in enumerate(sorted(FD,key=lambda x:-x[1])):D=Dl(sky(z,R,e,var));b.plot(l,D[2:31],color=plt.cm.viridis(j/3),lw=2,label=f'E={e:.3g} χ²={c2:.1f}, r²={r2(D):.2f}')
    for A in (a,b):A.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck');A.set_ylabel('D_ℓ (μK²)')
    a.set_yscale('log');a.set_ylim(1e-3,None);a.set_title(f'web riding CoE flow, {var}, R={R} Gpc, #{k}: '+(f'bracket guess E*≈{g["E"]:.3g}' if g['guess'] else f'no bracket (level-match E≈{ec:.3g})'))
    a.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85);b.set_ylim(0,2500);b.set_title(f'#{k} {var} R={R} Gpc: 4 closest level sets (linear)');b.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
  ax[-1].set_xlabel('ℓ');ax2[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig(O+f'flowR_{var}_ladder.png',dpi=90);fig2.tight_layout();fig2.savefig(O+f'flowR_{var}_lin4.png',dpi=90);plt.close('all')
