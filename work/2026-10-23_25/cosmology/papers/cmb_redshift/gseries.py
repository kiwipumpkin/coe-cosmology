"""Taylor coefficients of the universal LTB function g(w): R = M^(1/3) s^(2/3) g(w), w = E s^(2/3) / M^(2/3)
(G=c=1, s = t - t_B). Computed at 60 digits from the parametric solutions and a Chebyshev least-squares fit."""
import mpmath as mp, json
mp.mp.dps = 60
def g(w):
    w = mp.mpf(w)
    if w == 0: return mp.cbrt(mp.mpf(9)/2)
    T = (2*abs(w))**mp.mpf(1.5)
    if w > 0:
        eta = mp.findroot(lambda e: mp.sinh(e) - e - T, mp.cbrt(6*T))
        return (mp.cosh(eta) - 1)/(2*w)
    eta = mp.findroot(lambda e: e - mp.sin(e) - T, mp.cbrt(6*T))
    return (1 - mp.cos(eta))/(2*abs(w))
N = 26; a = mp.mpf('0.3')
xs = [a*mp.cos(mp.pi*(k+0.5)/80) for k in range(80)]
A = mp.matrix([[x**j for j in range(N)] for x in xs]); b = mp.matrix([g(x) for x in xs])
c = mp.lu_solve(A.T*A, A.T*b)
print([mp.nstr(c[j], 18) for j in range(8)])
test = [mp.mpf(v) for v in ('0.15', '-0.15', '0.05', '-0.2')]
for x in test: print(x, mp.nstr(sum(c[j]*x**j for j in range(N)) - g(x), 5))
json.dump([float(c[j]) for j in range(N)], open("gseries.json", "w"))
