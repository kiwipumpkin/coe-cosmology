import sys,os,json,numpy as np
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360));F=S.FLOOR
def mm(v):m=np.arange(1,13)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X=f['Xp'];M=f['Mp']
d=np.linalg.norm(X,axis=1);inn=d<S.RV;out=~inn
Pin=S.pattern(X[inn],M[inn]);Pout=S.pattern(X[out],M[out])
ring=S.RV*np.c_[np.cos(th),np.sin(th),0*th]
def shell_pattern(c):
    """exact uniform-shell (shell theorem) version of every outer body: mass M_j on a sphere of radius s_j=|x_j-c| about c.
    field at p = point mass at c if |p-c|>s_j, else 0. Same kernel convention as pattern(): [g(ring)-g(us)].n, mean removed."""
    s=np.linalg.norm(X[out]-c,axis=1);Mo=M[out]
    pr=np.linalg.norm(ring-c,axis=1);p0=np.linalg.norm(c)
    m_ring=(Mo[None,:]*(s[None,:]<pr[:,None])).sum(1)   # enclosed-in-sense: shells with s<|p-c| act as point at c
    m_us=Mo[s<p0].sum()
    if np.allclose(c,0): return np.zeros(360),m_ring,m_us
    unit=S.pattern(c[None,:],np.array([1.0]))  # can't use directly (mass differs per ring point) -> build by hand
    n=ring/S.RV;dv=c-ring;r3=np.linalg.norm(dv,axis=1)**3
    raw=m_ring*(dv*n).sum(1)/r3 - m_us*(c*n).sum(1)/np.linalg.norm(c)**3
    # normalise to pattern units using a single point mass check
    raw1=(dv*n).sum(1)/r3-(c*n).sum(1)/np.linalg.norm(c)**3; k=np.polyfit(raw1-raw1.mean(),unit,1)[0]
    v=k*raw;return v-v.mean(),m_ring,m_us
rows=[];out_json={}
for cn,c in (('us',np.zeros(3)),('CoE',np.asarray(S.C_G,float))):
    Pu,mr,mu=shell_pattern(c)
    out_json[cn]=dict(outer_uniform_only_m=mm(Pu).tolist(),outer_uniform_only_rms=float(S.rms(Pu)/F),mass_acting_ring_minmax=[float(mr.min()),float(mr.max())],mass_acting_us=float(mu))
    for fr in np.linspace(0,1,11):
        Po=(1-fr)*Pout+fr*Pu;res=S.R-Pin-Po;a=mm(res)
        rows.append(dict(center=cn,f_uniform=float(fr),m2=float(a[1]),total=float(S.rms(res)/F),per_m=a.tolist(),outer_m2=float(mm(Po)[1]),outer_total=float(S.rms(Po)/F)))
        print(cn,'%.1f'%fr,'m2 %.2f total %.2f | outer-only m2 %.2f rms %.2f'%(a[1],rows[-1]['total'],rows[-1]['outer_m2'],rows[-1]['outer_total']))
    print(cn,out_json[cn])
base=dict(no_web=mm(S.R).tolist(),no_web_total=float(S.rms(S.R)/F),inner_only=mm(S.R-Pin).tolist(),inner_only_total=float(S.rms(S.R-Pin)/F),inner_pattern_alone_m=mm(Pin).tolist(),outer_lumpy_alone_m=mm(Pout).tolist(),
 masses=dict(outer_net=float(M[out].sum()),inner_net=float(M[inn].sum()),outer_abs=float(abs(M[out]).sum()),n_in=int(inn.sum()),n_out=int(out.sum())))
print({k:(v[1] if isinstance(v,list) else v) for k,v in base.items() if k!='masses'},base['masses'])
json.dump(dict(base=base,shells=out_json,rows=rows),open('shell_exact.json','w'),indent=1)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(14,5.5))
for cn,c in (('us','C0'),('CoE','C3')):
    rr=[r for r in rows if r['center']==cn];x=[r['f_uniform'] for r in rr]
    ax[0].plot(x,[r['m2'] for r in rr],'o-',c=c,label=f'm=2, shell centred on {cn}');ax[0].plot(x,[r['outer_m2'] for r in rr],'s--',c=c,alpha=.6,label=f'outer-only m=2 ({cn})')
    ax[1].plot(x,[r['total'] for r in rr],'o-',c=c,label=f'total m=1-12 residual ({cn})');ax[1].plot(x,[r['outer_total'] for r in rr],'s--',c=c,alpha=.6,label=f'outer-only rms ({cn})')
for a,k,kt in ((ax[0],1,None),(ax[1],None,None)):
    pass
ax[0].axhline(base['inner_only'][1],c='k',ls=':',label='inner web only (outer removed)');ax[0].axhline(base['no_web'][1],c='gold',ls=':',label='no web (smooth ring)')
ax[1].axhline(base['inner_only_total'],c='k',ls=':');ax[1].axhline(base['no_web_total'],c='gold',ls=':')
for a in ax:a.set_xlabel('f = fraction of outer mass moved from lumpy web into uniform shells (0 lumpy, 1 shell-theorem limit)');a.set_yscale('log');a.legend(fontsize=8)
ax[0].set_ylabel('m=2 residual (× floor)');ax[1].set_ylabel('residual (× floor)')
fig.suptitle('Seed-3 web: outer (>14.3 Gpc) mass lumpy → uniform shells, inner 911 bodies fixed');plt.savefig('shell_sweep.png',dpi=110,bbox_inches='tight')
