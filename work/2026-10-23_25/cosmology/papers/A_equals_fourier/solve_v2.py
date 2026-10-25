"""Sedov-Taylor CoE analytic tidal model vs Migkas+2020 Fig 6 A(l): harmonic-amplitude solve,
tilt done properly (latitude-averaged, sector-averaged).  numpy only.
Assumptions: 90 deg longitude sectors (Migkas Fig 6: Delta l = 90, full latitude range), the
|b|<=20 deg Galactic-plane cut of the cluster sample, uniform sky weighting dOmega = cos b db dl.
mu = cos b cos b0 cos(l-l0) + sin b sin b0. The cut is symmetric in b, so sector averages factorise:
 <mu>   = a<c><C>,  <mu^2> = a^2<c^2><C^2> + s0^2<s^2>,  <mu^3> = a^3<c^3><C^3> + 3 a s0^2 <c s^2><C>
 (a=cos b0, s0=sin b0, c=cos b, s=sin b, C=cos(l-l0); <.>_b weighted by cos b over 20<|b|<90)."""
import numpy as np, json
d2r = np.pi/180
# ---------- latitude moments (closed form, weight cos b on [20,90] deg) ----------
bc = 20*d2r; hp = np.pi/2
F = lambda f: f(hp) - f(bc)
norm = F(np.sin)
c1 = F(lambda b: b/2 + np.sin(2*b)/4)/norm
c2 = F(lambda b: np.sin(b) - np.sin(b)**3/3)/norm
c3 = F(lambda b: 3*b/8 + np.sin(2*b)/4 + np.sin(4*b)/32)/norm
s2 = F(lambda b: np.sin(b)**3/3)/norm
cs2 = F(lambda b: b/8 - np.sin(4*b)/32)/norm
bb = np.linspace(bc, hp, 200001); ww = np.cos(bb)
num = lambda f: np.trapz(f*ww, bb)/np.trapz(ww, bb)
assert np.allclose([c1, c2, c3, s2, cs2], [num(np.cos(bb)), num(np.cos(bb)**2), num(np.cos(bb)**3),
                    num(np.sin(bb)**2), num(np.cos(bb)*np.sin(bb)**2)], atol=1e-8)
W = {m: np.sin(m*np.pi/4)/(m*np.pi/4) for m in (1,2,3)}
def coefs(b0deg):
    a, s0 = np.cos(b0deg*d2r), np.sin(b0deg*d2r)
    q0 = 1.5*a*a*c2 + 3*s0*s0*s2 - 1
    q2 = 1.5*a*a*c2*W[2]
    o1 = a*W[1]*(4.5*c1 - 45/8*a*a*c3 - 22.5*s0*s0*cs2)
    o3 = -15/8*a**3*c3*W[3]
    return dict(q0=q0, q2=q2, o1=o1, o3=o3)
# brute-force check of the factorised formulas (full 2-D sector average) at one b0
def brute(b0deg, N=1.0, ak=0.1, rho=0.5):
    b0 = b0deg*d2r
    bg = np.radians(np.arange(-89.875, 90, 0.25)); bg = bg[np.abs(bg) > bc]
    Lo = np.radians(np.arange(-44.75, 45, 0.5))
    B, L = np.meshgrid(bg, Lo, indexing='ij'); w = np.cos(B)
    psi = np.radians(np.arange(0, 360, 5.0)); out = []
    for p in psi:
        mu = np.cos(B)*np.cos(b0)*np.cos(p+L) + np.sin(B)*np.sin(b0)
        out.append(np.sum(w*N*(1+ak*(3*mu**2-1)+ak*rho*(4.5*mu-7.5*mu**3)))/np.sum(w))
    out = np.array(out)
    return out.mean(), [2*np.mean(out*np.cos(m*psi)) for m in (1,2,3)], [2*np.mean(out*np.sin(m*psi)) for m in (1,2,3)]
t = coefs(30); bf = brute(30)
print("check b0=30: analytic m0,m1,m2,m3 =", 1+0.1*t['q0'], 0.05*t['o1'], 0.1*t['q2'], 0.05*t['o3'], "\n   brute 2-D:", bf)
print(f"lat moments: <c>={c1:.5f} <c2>={c2:.5f} <c3>={c3:.5f} <s2>={s2:.5f} <cs2>={cs2:.5f}; W={W}")
# ---------- data ----------
dat = np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv', delimiter=',', names=True)
l_d, A_d = dat['l_deg'], dat['A_over_Aall']
X = np.column_stack([np.ones_like(l_d)] + [f(m*l_d*d2r) for m in (1,2,3) for f in (np.cos, np.sin)])
coef_refit = np.linalg.lstsq(X, A_d, rcond=None)[0]
coef = np.array([0.9965, 0.0296, 0.0434, 0.0636, 0.0366, -0.0299, -0.0114])
r2 = lambda pred: 1-np.sum((A_d-pred)**2)/np.sum((A_d-A_d.mean())**2)
R2_fit, R2_refit = r2(X@coef), r2(X@coef_refit)
amp = {0: coef[0]}; ph = {}
for m in (1,2,3):
    c, s = coef[2*m-1], coef[2*m]; amp[m] = np.hypot(c, s)
    ph[m] = (np.degrees(np.arctan2(s, c))/m) % (360/m)
print("N points", len(l_d), "refit", np.round(coef_refit, 4), "R2 given", round(R2_fit, 4), "refit", round(R2_refit, 4))
for m in (1,2,3): print(f"m{m}: amp={amp[m]:.5f}  m*phase={m*ph[m]:.3f}  peak l={ph[m]:.3f}")
l0 = ph[2]
# ---------- table ----------
b0s = np.arange(0, 91, 5.0)
print(" b0     q0       q2       o1       o3    |o1/o3|")
tab = []
for b in b0s:
    t = coefs(b); tab.append(dict(b0=b, **t))
    print(f"{b:4.0f} {t['q0']:8.4f} {t['q2']:8.4f} {t['o1']:8.4f} {t['o3']:8.4f} {abs(t['o1']/t['o3']) if abs(t['o3'])>1e-12 else np.inf:8.4f}")
# ---------- solve ----------
def bisect(f, a, b, tol=1e-12):
    fa = f(a)
    while b-a > tol:
        m = 0.5*(a+b); fm = f(m)
        if fa*fm <= 0: b = m
        else: a, fa = m, fm
    return 0.5*(a+b)
target = amp[1]/amp[3]
ratio = lambda b: abs(coefs(b)['o1']/coefs(b)['o3'])
fine = np.arange(0.05, 89.96, 0.05); rv = np.array([ratio(b) for b in fine])
roots = [bisect(lambda b: ratio(b)-target, fine[i], fine[i+1])
         for i in range(len(fine)-1) if (rv[i]-target)*(rv[i+1]-target) < 0]
imin = np.argmin(rv)
zero_o1 = [bisect(lambda b: coefs(b)['o1'], fine[i], fine[i+1]) for i in range(len(fine)-1)
           if coefs(fine[i])['o1']*coefs(fine[i+1])['o1'] < 0]
print(f"target |m1|/|m3|={target:.5f}; roots b0={roots}; min ratio {rv[imin]:.4f} at b0={fine[imin]:.2f}; o1=0 at b0={zero_o1}")
wrap = lambda x, P: (x + P/2) % P - P/2
L = l_d*d2r
sols = []
for b0 in roots:
    t = coefs(b0)
    Nak = amp[2]/abs(t['q2']); rho = amp[3]/(Nak*abs(t['o3'])); N = amp[0]-Nak*t['q0']; ak = Nak/N
    M = {0: N*(1+ak*t['q0']), 1: N*ak*rho*t['o1'], 2: N*ak*t['q2'], 3: N*ak*rho*t['o3']}
    nat = {m: (l0 + (0 if M[m] > 0 else 180/m)) % (360/m) for m in (1,2,3)}
    phi = {m: wrap(ph[m]-nat[m], 360/m) for m in (1,2,3)}
    def model(p1, p2, p3):
        return (M[0] + abs(M[1])*np.cos(L-(nat[1]+p1)*d2r) + abs(M[2])*np.cos(2*(L-(nat[2]+p2)*d2r))
                + abs(M[3])*np.cos(3*(L-(nat[3]+p3)*d2r)))
    s = dict(b0=b0, **t, Nak=Nak, rho=rho, N=N, ak=ak, Dr=1/rho, D=14.3/rho, Dgly=14.3/rho*3.2616,
             M0=M[0], M1=M[1], M2=M[2], M3=M[3], nat1=nat[1], nat2=nat[2], nat3=nat[3],
             phi1=phi[1], phi2=phi[2], phi3=phi[3], phi3_as_l3=3*phi[3],
             R2_nat=r2(model(0,0,0)), R2_rot=r2(model(phi[1],phi[2],0)), R2_all=r2(model(phi[1],phi[2],phi[3])))
    # Chris form: Amp*sin(l + phi):  |M1| cos(l - p) = |M1| sin(l + 90 - p)
    s['phiC_model'] = wrap(90-nat[1], 360); s['phiC_data'] = wrap(90-ph[1], 360)
    s['phi2C_model'] = wrap(90-2*nat[2], 360); s['phi2C_data'] = wrap(90-2*ph[2], 360)
    s['MoverMcrit'] = 2*ak   # if alpha = 1: kappa = M/(2 Mcrit)
    sols.append(s)
    for k, v in s.items(): print(f"  {k} = {v:.6f}")
psi = {}
for m in (1,2,3):
    c, sn = coef[2*m-1], coef[2*m]
    psi[m] = (c*np.cos(m*l0*d2r)+sn*np.sin(m*l0*d2r), -c*np.sin(m*l0*d2r)+sn*np.cos(m*l0*d2r))
    print(f"data m{m} psi frame C={psi[m][0]:+.5f} S={psi[m][1]:+.5f}  m*phase={np.degrees(np.arctan2(psi[m][1],psi[m][0])):.3f}")
print("ratio at b0=0,20,30,40,60:", [round(ratio(b),4) for b in (0.01,20,30,40,60)])
json.dump(dict(lat=dict(c1=c1,c2=c2,c3=c3,s2=s2,cs2=cs2), W=W, amp=amp, ph=ph, l0=l0, R2_fit=R2_fit,
               R2_refit=R2_refit, coef_refit=list(coef_refit), target=target, roots=roots, table=tab, sols=sols,
               psi=psi), open('solve_v2_results.json', 'w'), indent=1, default=float)
