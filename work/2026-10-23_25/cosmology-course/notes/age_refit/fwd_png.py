import json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,numpy as np,textwrap
J=json.load(open('fwd/jerk.json'));R=json.load(open('fwd/eta_vc_rrec.json'))
def g(l,s=''): return json.load(open(f'fwd/{l}{s}.json'))['free']
p1=["FORWARD REFIT (forward.py) — page 1: fits",
"Free: M [kg], G_i=2GM_enc/(r_i V0^2), p (V0~r0^p), f=matter/M inside r_i, L, shape (n/alpha/eta), l0.",
"Derived: V0, D_today, H_par, H_perp, age, q, j. Start r_i fixed (exact scale symmetry), then rescaled so mean H=70.",
f"JERK (local Taylor, H_par/H_perp,k2,k3,D,l0): chi2 {J['chi2']:.2f}, R2 {J['R2']:.2f}, {J['sigma_k5']:.1f} sigma, D {J['D_Gpc']:.2f} Gpc, H_par {J['Hpar']:.0f}, H_perp {J['Hper']:.0f}"]
for s,lab in (('_vc','v_o <= c'),('','no speed limit')):
    p1.append(f"--- {lab} ---")
    for l in ('exp','soft','eta','pl'):
        d=g(l,s); p1.append(f"{l}: chi2 {d['chi2']:.2f} R2 {d['R2']:.2f} {d['sigma']:.1f}sig | M {d['M_kg']:.1e} kg | D {d['D_today_Gpc']:.3g} Gpc | v/c {d['vo_over_c']:.3g} | H {d['H_par']:.1f}/{d['H_perp']:.1f} | age {d['age_Gyr']:.1f} | grav/kin {d['grav_over_kin_today']:.3f}")
e=g('eta','_vc'); r=R['rows'][1]
p2=["page 2: CMB checks (same universe, eta v<=c best)",
f"Temperature-traced LS (T=3000 K): found {e['cmb_found']}/360, mean T_obs {e['Tobs_mean']:.0f} K (needs 2.725), cancel {e['cancel_pct']:.0f}%, resid {e['resid_over_1e5']:.0f}x floor",
f"t_LS: 0deg {e['tLS_0']:.2f}, 180deg {e['tLS_180']:.2f} Gyr; shell crossing = {e['shell_crossing']}",
"Chris R_rec picture: plasma = |x|<R_rec; LS on near side of sphere, only within psi_crit=asin(R_rec/D) of CoE.",
"gap(180) = integral dr/U - (D-R_rec)/c  ~  (D-R_rec)(1/U_o - 1/c): faster U_o closes the gap",
"Closed form 0 vs 180: CoE frame dt=2D/c;  observer frame dt'=(2/c)[D sqrt(1+U^2/c^2) - U R_rec/c]"]
for r in R['rows']: p2.append(f"R_rec/D {r['Rrec_over_D']}: psi_crit {r['psi_crit_deg']:.0f}deg, R_crit {r['R_critical_Gpc']:.3f} Gpc, gap int {r['gap180_int']:.2f} vs constU {r['gap180_constU']:.2f} Gyr, dt' {r['dt_formula_obs']:.2f}")
for i,p in enumerate((p1,p2),1):
    fig=plt.figure(figsize=(16,10));y=0.97
    for t in p:
        for w in textwrap.wrap(t,95): fig.text(0.02,y,w,fontsize=15,family='monospace',va='top');y-=0.037
    fig.savefig(f'fwd/forward_summary_p{i}.png',dpi=90);plt.close()
