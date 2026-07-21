"""
NOTE 13, stage A-probe (self-contained): WHO IS LYING -- symbolic A_L or trackers?
1) recompute symbolic layer coefficients AL, BL (odd d, eps = 1/d), print certs.
2) dense scan -> ALL real roots of p_d(t) = s in the pin window, 100-digit polish,
   all tau; show which root binds (min tau) and which family AL/BL describes.
"""
import json
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 100
ORD = 3
x, Y, s = sp.symbols("x Y s")
S5 = (s + 5) / 2
L0 = sp.log(S5)

def pow_series(shift, order):
    lser = sp.series(sp.log(1 + x*Y), Y, 0, order).removeO().expand()/Y
    e    = sp.series(sp.exp(lser), Y, 0, order).removeO().expand()
    fac  = sp.series((1 + x*Y)**shift, Y, 0, order).removeO().expand()
    return sp.expand(e * fac)

E1   = pow_series(0, ORD)
Em2  = sp.expand(pow_series(0, ORD) * sp.series((1+x*Y)**(-2), Y, 0, ORD).removeO())
E1p1 = sp.expand(pow_series(1, ORD))

t   = -1 - x*Y
Yo1 = sp.series(Y/(1+Y), Y, 0, ORD+1).removeO().expand()
cd  = sp.expand(6*Y*Yo1)
p   = 2*t - 3*t**2 + t*(1 - t)*((-Em2) - cd)
Phi = t**2 - t**3 + Y*(-1)*E1 - E1p1*Yo1 - cd*(t**2/2 - t**3/3)

x0, x1, x2 = sp.symbols("x0 x1 x2")
xsub = x0 + x1*Y + x2*Y**2
E0 = sp.symbols("E0")

def Yseries(e):
    e = sp.expand(e.subs(x, xsub))
    for a in list(e.atoms(sp.exp)):
        arg  = sp.expand(a.args[0]); arg0 = arg.subs(Y, 0)
        rest = sp.expand(arg - arg0)
        ser  = sp.series(sp.exp(rest), Y, 0, ORD).removeO().expand()
        e = sp.expand(e.subs(a, (E0 if arg0 == x0 else sp.exp(arg0))*ser))
    return sp.expand(e)

def cY(e, k): return sp.expand(e).coeff(Y, k)

pexp = Yseries(p - s)
eq1, eq2 = cY(pexp, 1), cY(pexp, 2)
sol1 = sp.solve(sp.Eq(eq1, 0), x1)[0]
sol2 = sp.solve(sp.Eq(eq2.subs(x1, sol1), 0), x2)[0]
tauexp = Yseries(t*s - Phi)
# x1, x2 are SOLVED perturbations: substitute them into the tau series, THEN
# mark the pin constants E0 -> (s+5)/2, x0 -> ln((s+5)/2).
tau0 = sp.simplify(cY(tauexp, 0).subs(E0, S5).subs(x0, L0))
tau1 = sp.simplify(cY(tauexp, 1).subs(x1, sol1).subs(E0, S5).subs(x0, L0))
tau2 = sp.simplify(cY(tauexp, 2).subs([(x1, sol1), (x2, sol2)]).subs(E0, S5).subs(x0, L0))
ALc  = (s + 5) - 5*sp.log((s+5)/2)
print("CERT tau0 = -s-2 exact:", sp.simplify(tau0 + s + 2) == 0)
print("CERT tau1 = (s+5)-5 ln((s+5)/2) exact:", sp.simplify(tau1 - ALc) == 0)
print("tau1 =", tau1)
print("tau2 =", tau2)
json.dump({"AL": str(tau1), "BL": str(tau2),
           "cert0": bool(sp.simplify(tau0 + s + 2) == 0),
           "cert1": bool(sp.simplify(tau1 - ALc) == 0)},
          open("ring_probe_sym.json", "w"), indent=1)
ALf = sp.lambdify(s, tau1, "math"); BLf = sp.lambdify(s, tau2, "math")

# ---------------- dense scans ----------------
def seed_co(d):
    w = sp.symbols("w")
    p_ = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - sp.Rational(6, d*(d+1))))
    Ph = sp.expand(sp.integrate(p_, w)); ta = sp.expand(w*p_ - Ph)
    return ([mp.mpf(str(c)) for c in sp.Poly(p_, w).all_coeffs()],
            [mp.mpf(str(c)) for c in sp.Poly(ta, w).all_coeffs()])

def all_roots(sv, pc, tc):
    pcs = list(pc); pcs[-1] -= mp.mpf(sv)
    dps = [c*(len(pcs)-1-i) for i, c in enumerate(pcs[:-1])]
    f, df = (lambda z: mp.polyval(pcs, z)), (lambda z: mp.polyval(dps, z))
    grid = np.concatenate([np.linspace(-1.7, -0.9, 6001), np.linspace(-0.9, 2.5, 2501)])
    outs = []
    for i in range(len(grid)-1):
        a, b = mp.mpf(float(grid[i])), mp.mpf(float(grid[i+1]))
        if f(a)*f(b) < 0:
            try:
                z = mp.findroot(f, (a, b), solver="bisect", tol=mp.mpf(10)**(-60))
            except ValueError:
                continue
            if all(abs(z - q) > 1e-25 for q, _ in outs):
                outs.append((mp.re(z), mp.polyval(tc, z)))
    return sorted(outs, key=lambda q: q[1])

for d in (61, 201):
    pc, tc = seed_co(d)
    for sv in (-0.5, 0.3, 0.0, -1.2, -2.0):
        rr = all_roots(sv, pc, tc)
        a, b = ALf(sv), BLf(sv)
        print(f"\n d={d} s={sv:+.1f}:  {len(rr)} real roots;  A_L={a:.8f} B_L={b:+.6f}  A_L+B_L/d={a+b/d:.8f}")
        for z, tau in rr:
            print(f"   t = {mp.nstr(z, 18)}   tau = {mp.nstr(tau, 18)}   d*(tau+s+2) = {float(d*(tau+mp.mpf(sv)+2)):+.8f}")
print("\nsaved ring_probe_sym.json")
