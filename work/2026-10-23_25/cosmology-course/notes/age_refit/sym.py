"""Odd/even decomposition of visible dT/T about the CoE axis (model angle th: 0=away from CoE, 180=CoE)."""
import json,glob,sys,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
def load(f):
    d=json.load(open(f)); R=d['rows']
    if 'z_fwd' in str(R[:3]) or any('z_fwd' in r for r in R): pts=[(r['ell'],r['z_fwd']) for r in R if r.get('z_fwd') is not None]; kind='common-event'
    else: pts=[(r['ell'],r['z']) for r in R if r.get('status')=='rec']; kind='3000K'
    return d,kind,pts
res=[]
for f in sorted(glob.glob('gr2/*_*_*.json'))+sorted(glob.glob('gr2/rec3000_*.json')):
    if 'validation' in f or f.startswith('gr2/rec3000')!=('rec3000' in f): continue
    d,kind,pts=load(f)
    if len(pts)<6: res.append(dict(file=f,kind=kind,visible=len(pts),note='too few visible directions')); continue
    th=np.array([p[0] for p in pts]); S=1/(1+np.array([p[1] for p in pts])); D=S/S.mean()-1
    A=np.c_[np.ones(len(th)),np.cos(np.radians(th)),np.sin(np.radians(th))]; c=np.linalg.lstsq(A,np.log(S),rcond=None)[0]
    Sn=S/np.exp(A[:,1:]@c[1:]); Dn=Sn/Sn.mean()-1
    out=dict(file=f,tag=d['tag'],kind=kind,visible=len(pts),total=len(d['rows']),m1_fit_cos=float(c[1]),m1_fit_sin=float(c[2]))
    for lab,v in (('raw',D),('m1removed',Dn)):
        m=dict(zip(np.round(th,6),v)); pr=[(t,m[t],m.get(round((360-t)%360,6))) for t in th if t<=180 and round((360-t)%360,6) in m]
        a=np.array([p[1] for p in pr]); b=np.array([p[2] for p in pr]); odd=(a-b)/2; even=(a+b)/2
        left=[t for t in th if t<180]; right=[t for t in th if t>180]
        tl,tr_=max(left),min(right)
        out[lab]=dict(rms_total=float(np.sqrt(np.mean(v**2))),rms_odd=float(np.sqrt(np.mean(odd**2))),rms_even=float(np.sqrt(np.mean(even**2))),max_abs_odd=float(np.abs(odd).max()),
            p2p_even=float(np.ptp(even)),edge_left=(float(tl),float(m[tl])),edge_right=(float(tr_),float(m[tr_])),jump_at_180=float(m[tr_]-m[tl]),
            gap_deg=float(tr_-tl), value_at_180=float(m[180.0]) if 180.0 in m else None)
        out[lab+'_curves']=(list(map(float,[p[0] for p in pr])),list(map(float,a)),list(map(float,b)),list(map(float,even)))
    res.append(out)
    fig,ax=plt.subplots(2,1,figsize=(16,11))
    for k,lab in enumerate(('raw','m1removed')):
        t,a,b,e=out[lab+'_curves']
        ax[k].plot(t,a,'o-',lw=3,ms=8,label='0°→180° half f(θ)'); ax[k].plot(t,b,'s--',lw=3,ms=8,label='360°→180° half f(360°−θ)'); ax[k].plot(t,e,'k-',lw=2,label='even part (average of halves)')
        ax[k].plot(t,(np.array(a)-np.array(b))/2,'r:',lw=3,label='odd part')
        o=out[lab]; ax[k].set_title(f"{lab}: rms odd {o['rms_odd']:.2e}, rms even {o['rms_even']:.2e}, jump at 180° {o['jump_at_180']:.2e} (gap {o['gap_deg']:.0f}°)",fontsize=16)
        ax[k].set_ylabel('dT/T',fontsize=18); ax[k].legend(fontsize=14); ax[k].tick_params(labelsize=14); ax[k].grid(); ax[k].set_xlim(0,180)
    ax[1].set_xlabel('angle from anti-CoE direction θ (deg); 180 = CoE',fontsize=18)
    fig.suptitle(f"{kind} {d['tag']}: visible {len(pts)}/{len(d['rows'])}",fontsize=18); fig.tight_layout()
    fig.savefig(f"gr2/sym_{'3000_' if kind=='3000K' else ''}{d['tag']}.png",dpi=80); plt.close(fig)
for r in res: [r.pop(k) for k in list(r) if k.endswith('_curves')]
json.dump(res,open('gr2/symmetry.json','w'),indent=1)
for r in res:
    if 'raw' not in r: print(r); continue
    print(f"{r['kind']:12s} {r['tag']:14s} vis {r['visible']}/{r['total']} | raw: odd {r['raw']['rms_odd']:.1e} even {r['raw']['rms_even']:.3g} jump180 {r['raw']['jump_at_180']:.3g} edges {r['raw']['edge_left']},{r['raw']['edge_right']} | m1rem: odd {r['m1removed']['rms_odd']:.1e} even {r['m1removed']['rms_even']:.3g} p2p {r['m1removed']['p2p_even']:.3g} jump {r['m1removed']['jump_at_180']:.3g}")
