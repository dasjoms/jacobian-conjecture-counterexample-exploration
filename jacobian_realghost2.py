"""
Note 5, stage 2: design a seed p (deg 4, rational coeffs) whose Keller lift is
surjective over C (no (3,2) / (5) all-multiple patterns) AND over R
(every real cusp of Phi has real simple partners).

Recipe recap (c = b = 1): gamma = 1 + a*v + t, u = 1+v, w = u*gamma,
F = (u + q(w)/g^2)/x^2, (1 + p(w)/g)/x, x*g , q = w*p - int(p), a = -(1+kap)/(2+kap).
Fiber equation (C != 0): Phi(w) = BC*w - AC^2 ; escape <-> multiple root.
"""
import sympy as sp
from sympy import symbols, Rational as R, diff, integrate, expand, factor, cancel
import itertools, random

w = symbols("w")

def outer_quadratics(c4, c3, c2, c1):
    """For each real critical point r of p, divide h = Phi - s0*w + r0 by (w-r)^3
    by exact high-precision long division; return list of (r, quotient coefficients,
    division remainder norm)."""
    p = c4*w**4 + c3*w**3 + c2*w**2 + c1*w
    Phi = expand(integrate(p, w))
    pp = diff(p, w)
    crits = sp.nroots(pp, n=40, maxsteps=200)
    crits = [r for r in crits if abs(sp.im(r)) < sp.Float('1e-30')]
    out = []
    for r in crits:
        s0 = p.subs(w, r); r0 = r*s0 - Phi.subs(w, r)
        h = Phi - s0*w + r0       # exact rational in coefficients of Float(r)
        a = [sp.Float(h.expand().coeff(w, k), 50) for k in range(5, -1, -1)]
        r = sp.Float(r, 50)
        d = [sp.Float(1, 50), -3*r, 3*r**2, -r**3]
        # long division of a (len 6) by d (len 4) -> quotient len 3, remainder len 3
        rem = a[:]
        quo = []
        for k in range(3):
            qk = rem[k]
            quo.append(qk)
            for i in range(4):
                rem[k+i] = rem[k+i] - qk*d[i]
        resid = rem[3:]
        out.append((float(r), quo, [complex(v) for v in resid]))
    return out

def cusp_check(c4, c3, c2, c1, verbose=False):
    """ok=True iff every REAL critical point of p has outer quadratic
    with positive discriminant (2 real simple partners)."""
    ok = True
    details = []
    try:
        rows = outer_quadratics(c4, c3, c2, c1)
    except Exception as e:
        return False, [("nroots failed", str(e), None)]
    for r, quo, resid in rows:
        q2, q1, q0 = quo
        disc = q1**2 - 4*q2*q0
        resid_norm = max(abs(v) for v in resid)
        this_ok = disc > 0 and resid_norm < 1e-20
        details.append((r, float(disc), resid_norm, this_ok))
        ok = ok and this_ok
    return ok, details

random.seed(42)
found = []
tried = 0
for c4, c3, c2, c1 in itertools.product(range(-3,4), range(-3,4), range(-3,4), range(-3,4)):
    if c4 == 0: continue
    p = c4*w**4 + c3*w**3 + c2*w**2 + c1*w
    kap = diff(p, w).subs(w, 1)
    if kap <= -3 or kap == -2 or c1 == 0:   # want gamma-sheet structure too (p'(0)!=0); skip kappa near -2
        continue
    tried += 1
    # quick pre-filter: p' must have exactly 1 real root region check skipped; do full check
    ok, det = cusp_check(c4, c3, c2, c1)
    if ok and det:   # has at least one real cusp, all rescued
        found.append((c4, c3, c2, c1, kap, det))
    if len(found) >= 8:
        break
print(f"tried {tried} candidates; {len(found)} pass all real cusps:")
for c4, c3, c2, c1, kap, det in found:
    print(f"  p = {c4}w^4 + {c3}w^3 + {c2}w^2 + {c1}w   kappa = {kap}")
    for r0, ok, outer in det:
        print(f"    cusp at w0={r0:.6f}: outer roots {outer} rescued={ok}")
