"""
Note 5, stage 2b: search the CORRECT design space.
Normalization for the lift to be polynomial (b = c = 1):
    p(1) = -1,  Phi(1) = int_0^1 p = 0,  kappa = p'(1) != -2.
Free parameters: (c2, c3, c4) of p = c0 + c1 w + c2 w^2 + c3 w^3 + c4 w^4 with
    c1 = 2*(S1 - S2); c0 = 2*S2 - S1;
    S1 = -1 - c2 - c3 - c4;  S2 = -(c2/3 + c3/4 + c4/5).
Want: every real cusp of Phi (root of p') rescued, i.e. the outer quadratic
of h = Phi - p(w0) w + (w0 p(w0) - Phi(w0)) has positive discriminant.

Discriminant closed form (from matching (w-w0)^3 (a w^2 + b w + g) to h):
    a = c4/5, b = c3/4 + 3 w0 c4/5, g = c2/3 + 3 w0 c3/4 + 6 w0^2 c4/5
    D = c3^2/16 - (3/10) c3 c4 w0 - (3/5) c4^2 w0^2 - (4/15) c4 c2.
"""
import sympy as sp
from sympy import symbols, Rational as R, diff, integrate, expand
import itertools

w = symbols("w")

def seed(c2, c3, c4):
    S1 = -1 - c2 - c3 - c4
    S2 = -(R(c2,1)/3 + R(c3,1)/4 + R(c4,1)/5)
    c1 = 2*(S1 - S2); c0 = 2*S2 - S1
    p = [sp.nsimplify(c0), sp.nsimplify(c1), sp.nsimplify(c2), sp.nsimplify(c3), sp.nsimplify(c4)]
    return p

def cusp_data(pco):
    c0,c1,c2,c3,c4 = pco
    p = sum(ci*w**i for i,ci in enumerate(pco))
    Phi = expand(integrate(p, w))
    pp = diff(p, w)
    crits = sp.nroots(pp, n=40, maxsteps=200)
    crits = [r for r in crits if abs(sp.im(r)) < sp.Float('1e-28')]
    rows = []
    for r in crits:
        r50 = sp.Float(r, 50)
        D = (R(c3)**2/16 - R(3,10)*c3*c4*r50 - R(3,5)*R(c4)**2*r50**2 - R(4,15)*c4*c2)
        rows.append((float(r50), sp.simplify(D)))
    return rows

hits = []
count = 0
for c2, c3, c4 in itertools.product(range(-4,5), range(-4,5), range(-4,5)):
    if c4 == 0: continue
    pco = seed(c2, c3, c4)
    p = sum(ci*w**i for i,ci in enumerate(pco))
    kap = sp.diff(p, w).subs(w, 1)
    if kap == -2: continue
    if abs(kap + 2) < sp.Rational(1,2): continue
    count += 1
    try:
        rows = cusp_data(pco)
    except Exception:
        continue
    if rows and all(sp.N(D, 20) > 0 for _, D in rows):
        hits.append((pco, kap, rows))
print(f"searched {count} seeds; {len(hits)} fully rescued")
for pco, kap, rows in hits[:12]:
    print("  p(w) =", sum(ci*w**i for i,ci in enumerate(pco)).expand(), " kappa =", kap,
          " cusp D's:", [(round(r,4), sp.nsimplify(D)) for _, D in rows])
