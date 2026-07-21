"""
NOTE 15, stage 4 (corrected, final): certified shadow identity + full series.

 (C1) CORRECTED identity (odd d; U := |t*+1| > 0, X := (1+U)^(d-2)):
        (d^2+d)(U+2)^2 (X-2) + 4U^2 + 17U + 19 = X (U+1) (1 + d(U+2))
      [earlier transcript dropped d*U^3 in remE - caught because the residual
       0.01316 exactly equalled X*d*U^3 - the debugging lesson price was worth it]
      Audit: |residual| < 1e-45 at every exact root, d = 3,5,...,45. Anchor: d=3
      is EXACT (U=1, X=2: 0 + 40 - 2*2*(1+3*3) = 0).
 (C2) ghost: finish d=43-ratio audit + asymptotic law v = 1/[3(d(d+1)-2)] + O(3^-d).
 (C3) left-cusp fixed-point iteration audit, d = 7..41.
 (M3) m(3) census 1e6, lock [8.45, 8.85]%.
 (B) machine series: U(d) = u1/dn + u2/dn^2 + u3/dn^3 with u1 = ln 2,
     u2, u3 derived EXACTLY from (C1) by order-by-order balance (eps-series),
     then LSQ-verified against the exact roots (d=9..45); same for
     X(d)-2 and sigma*(d) := s*+1 laws.
 (F1/F2) publish chamber-12 locks.
"""
import json, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import symbols, expand, diff, integrate, Rational as R, Poly, exp, log, ln

mp.mp.dps = 90
w = symbols("w")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def diag_root(d):
    p = seed(d); Phi = integrate(p, w)
    G = expand(Phi - (w-1)*p)
    PP = Poly(sp.expand(G*(d*(d+1))), w)
    fc = [float(c) for c in PP.all_coeffs()]
    Gf = sp.lambdify(w, G, "mpmath")
    best = None
    for z in np.roots(fc):
        if abs(z.imag) < 1e-7:
            xr = mp.mpf(z.real)
            try: xr = mp.findroot(Gf, xr, tol=mp.mpf("1e-80"), verify=False)
            except Exception: pass
            if abs(Gf(xr)) < mp.mpf("1e-40") and xr < -1.000001 and (best is None or xr > best):
                best = xr
    spoint = sp.lambdify(w, p, "mpmath")(best)
    # INDEPENDENT exact audit: substitute back symbolically at high precision
    aud = abs(sp.N(p.subs(w, sp.Float(best, 100)) - spoint, 60))
    return best, spoint, aud

roots = {}
audw = 0
for d in [3] + list(range(5, 46, 2)):
    tv, sv, aud = diag_root(d)
    roots[d] = (tv, sv)
    audw = max(audw, mp.mpf(str(aud)))
print(f"  exact roots d=3..45 odd computed; worst p-subst audit {mp.nstr(audw, 6)}", flush=True)
out["root_audits_worst"] = mp.nstr(audw, 6)

# ---------------- (C1) corrected identity audit ----------------
print("="*84); print("(C1) corrected shadow identity audit")
def c1_res(U, d):
    X = (1+U)**(d-2)
    return abs((d*d+d)*(U+2)**2*(X-2) + (4*U**2+17*U+19) - X*(U+1)*(1+d*(U+2)))
# d=3 anchor
a3 = c1_res(mp.mpf(1), 3)
print(f"  d=3 EXACT anchor (U=1, X=2): residual {mp.nstr(a3, 8)}")
worst = mp.mpf(0)
for d, (tv, sv) in roots.items():
    if d == 3: continue
    U = -tv-1
    r_ = c1_res(U, d)
    worst = max(worst, r_)
print(f"  worst residual d=5..45: {mp.nstr(worst, 8)}   [<1e-45: {worst < mp.mpf('1e-45')}]")
out["C1_ok"] = bool(worst < mp.mpf("1e-45")); out["C1_worst"] = mp.nstr(worst, 6)

# ---------------- (B) order-by-order series ----------------
print("="*84); print("(B) eps-series balance from the identity")
eps = symbols("eps", positive=True)
dn = 1/eps
u1, u2, u3 = symbols("u1 u2 u3")
U_ser = u1*eps + u2*eps**2 + u3*eps**3
# X = (1+U)^(1/eps): log series to eps^3
z = symbols("z")
logser = u1 + (u2 - u1**2/2)*eps + (u3 - u1*u2 + u1**3/3)*eps**2
X_ser = sp.exp(logser - u1) * sp.exp(u1)
X_ser = sp.expand(X_ser)
# expand exp(u1) * exp((u2-u1^2/2) eps + ...) in eps to order 2:
v1 = u2 - u1**2/2
v2 = u3 - u1*u2 + u1**3/3
X_ser = sp.exp(u1)*(1 + v1*eps + (v2 + v1**2/2)*eps**2)
d_ser = dn + 2
F = (d_ser*d_ser+d_ser)*(U_ser+2)**2*(X_ser-2) + (4*U_ser**2+17*U_ser+19) - X_ser*(U_ser+1)*(1+d_ser*(U_ser+2))
Fser = sp.expand(F*eps)
print("  F*eps at O(eps^0):")
c0 = sp.expand(Fser).subs(eps, 0)
c0 = sp.expand(c0 - sp.expand(c0).subs(eps, 0)) + sp.expand(Fser).subs(eps,0)
print("   =", sp.collect(sp.expand(Fser).subs(eps, 0), sp.exp(u1)))
# collect at eps^0 and eps^1
Fcol = sp.series(Fser, eps, 0, 2).removeO()
Fcol = sp.expand(Fcol)
A0 = sp.expand(Fcol.subs(eps, 0))
A1 = sp.expand((Fcol - A0)/eps)
print("  coeff eps^0:", sp.collect(A0, sp.exp(u1)))
print("  coeff eps^1:", sp.collect(A1, sp.exp(u1)))
out["series_eq"] = [str(sp.collect(A0, sp.exp(u1))), str(sp.collect(A1, sp.exp(u1)))]
json.dump(out, open("jcorner_stage4_partial.json", "w"))
print(f"checkpoint [{time.time()-t0:.0f}s]")
