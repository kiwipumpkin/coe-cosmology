import json,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,textwrap
def J(p):
    try: return json.load(open(p))
    except Exception: return None
L=[];F='%.4g'
L.append("SIGN-FIX + GR RERUN: before (old sign Phi_o-Phi_e, unoriented global CMB) | Newton fixed | exact Schwarzschild (Cosmology)")
L.append("Migkas fits: no potential term -> identical chi2 before/after (rerun reproduced to all digits).")
L.append("law/limit: chi2 | rs/D (enc) | rs/D (pt) | global-traced cancel%  before / newton / GR | Tmean K before/newton/GR")
for s,lab in (('_vc','v<=c'),('','nolim')):
    for l in ('exp','soft','eta','pl'):
        b,n,g=[(J(f'fwd/{d}/{l}{s}.json') or {}).get('free') for d in ('before','after_newton','after_gr')]
        if not n: continue
        rsD=n['grav_over_kin_today']*n['vo_over_c']**2; rsp=n['grav_over_kin_pointmass_today']*n['vo_over_c']**2
        flag=' INSIDE r_s!' if rsD>1 else ''
        c=lambda x,k: (F%x[k]) if x and k in x and x[k] is not None else '-'
        L.append(f"{l} {lab}: chi2 {n['chi2']:.2f} | {rsD:.2e} | {rsp:.2e}{flag} | {c(b,'cancel_pct')} / {c(n,'cancel_pct')} / {c(g,'cancel_pct')} | {c(b,'Tobs_mean')}/{c(n,'Tobs_mean')}/{c(g,'Tobs_mean')}")
L.append("LOCAL-SR CMB (q=R/D 0.23; cancel% phi=l0 / phi free): before | newton | GR")
for l in ('exp','soft','eta','pl'):
    r=[(J(f'fwd/{d}/sr_{l}.json') or {}).get('_vc',{}).get('0.23',{}) for d in ('before','after_newton','after_gr')]
    L.append(f"{l} v<=c: "+" | ".join((f"{x['cancel_l0']:.2f} / {x['cancel_phifree']:.2f}" if x.get('ok') else f"no full sky ({x.get('nfound')}/360)") for x in r))
for nm,key in (('bracket','cancel'),('bracket_sr','cancel_phifree')):
    for d in ('before','after_newton','after_gr'):
        x=J(f'fwd/{d}/{nm}.json')
        if not x: continue
        rows=[r for r in x['rows'] if r.get('ok') and np.isfinite(r.get(key) if r.get(key) is not None else np.nan)]
        if not rows: L.append(f"{nm} {d}: no valid rows"); continue
        B=rows[-1]; best=max(rows,key=lambda r:r[key])
        L.append(f"{nm} {d}: B end chi2 {B['chi2']:.1f} cancel {B[key]:.1f}% M {B['M_kg']:.2g} kg gk {B['gk']:.3g} | best cancel on path {best[key]:.1f}% at f={best['f']:.2f} (chi2 {best['chi2']:.1f})")
s2=J('fwd/step2_redo.json')
if s2:
    L.append("STEP 2 same-universe (unified2): cancel% before -> after sign fix")
    L.append(" ; ".join(f"{k} {v['before_cancel']:.1f}->{v['cancel']:.1f}" for k,v in s2.items()))
open('fwd/compare_fix.txt','w').write("\n".join(L)); print("\n".join(L))
lines=[w for t in L for w in textwrap.wrap(t,92)]; per=24
for p in range(0,len(lines),per):
    fig=plt.figure(figsize=(16,10)); y=0.97
    for w in lines[p:p+per]: fig.text(0.015,y,w,fontsize=15,family='monospace',va='top'); y-=0.04
    fig.savefig(f'fwd/compare_fix_p{p//per+1}.png',dpi=90); plt.close()
