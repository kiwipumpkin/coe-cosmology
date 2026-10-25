import json,numpy as np,textwrap,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
L=[]
def g(d,*k):
    for kk in k:
        if not isinstance(d,dict) or kk not in d: return np.nan
        d=d[kk]
    return d
L.append("LIGHT-CONE MIGKAS + BACKWARD + UNIFORM-SURFACE CMB (GR = exact Schwarzschild; Newton in [])")
L.append("law mode: chi2 (old today-pos same params) R2 sig | M kg | D Gpc | v/c | H||/Hperp | age | grav/kin | rs/D")
for m in ('fwd','back'):
  for l in ('exp','soft','eta','pl'):
    d=json.load(open(f'fwd/lc_{l}_{m}_gr.json'))
    ex=f" origin={d['reached_origin']} bang-spread {d['bang_time_spread_Gyr']:.2f}Gyr" if m=='back' else ''
    L.append(f"{l} {m}: {d['chi2']:.2f} ({d['chi2_old_todaypositions_same_params']:.1f}) {d['R2']:.2f} {d['sigma']:.1f}s | {d['M_kg']:.2g} | {d['D_Gpc']:.3g} | {d['vo_c']:.3f} | {d['Hpar']:.0f}/{d['Hper']:.0f} | {d['age_Gyr']:.1f} | {d['gk']:.3f} | {d['rs_over_D']:.1e}{ex}")
L.append("CMB cancel% vs Migkas curve (phi=l0 / phi free), raw | m=1 removed.  GR [Newton]")
for m in ('fwd','back'):
  for l in ('exp','soft','eta','pl'):
    d=json.load(open(f'fwd/lc_{l}_{m}_gr.json')); n=json.load(open(f'fwd/lc_{l}_{m}.json'))
    c=d['cmb'].get('0.23',{}); cn=n['cmb'].get('0.23',{}); u=d['cmb_uniform']['1100']; un=n['cmb_uniform']['1100']
    f=lambda x,k: f"{x[k]:.1f}" if k in x else 'n/a'
    L.append(f"{l} {m} localSR q.23: {f(c,'cancel_l0')}/{f(c,'cancel_free')} | {f(c,'cancel_l0_nodip')}/{f(c,'cancel_free_nodip')} [{f(cn,'cancel_l0')}|{f(cn,'cancel_l0_nodip')}]")
    L.append(f"   uniform a=1/1100 (z~{u['mean_z']:.3g}, inside {u['n_inside']}/360): {f(u,'cancel_l0')}/{f(u,'cancel_free')} | {f(u,'cancel_l0_nodip')}/{f(u,'cancel_free_nodip')} [{f(un,'cancel_l0')}|{f(un,'cancel_l0_nodip')}]  dipole {g(u,'S_modes','1'):.3g} quad {g(u,'S_modes','2'):.3g}")
open('fwd/lc_summary.txt','w').write("\n".join(L)); print("\n".join(L))
ls=[w for t in L for w in textwrap.wrap(t,100)]; per=24
for p in range(0,len(ls),per):
    fig=plt.figure(figsize=(16,10)); y=0.97
    for w in ls[p:p+per]: fig.text(0.012,y,w,fontsize=14,family='monospace',va='top'); y-=0.04
    fig.savefig(f'fwd/lc_summary_p{p//per+1}.png',dpi=90); plt.close()
