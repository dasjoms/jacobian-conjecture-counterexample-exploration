"""
Stage 2c: WHY do all cusps come out unrescued (D < 0) in normalized space?
Look at the joint distribution of (#real cusps, signs of D) over the grid,
plus a symbolic look at D as a quadratic form under the normalization.

Recall: p = c0 + c1 w + c2 w^2 + c3 w^3 + c4 w^4 (deg 4), normalized by
  p(1) = -1, Phi(1) = 0   (lift polynomiality), so 3 free params (c2,c3,c4).
At a cusp w0 (real root of p'), the outer quadratic has discriminant
  D = c3^2/16 - (3/10) c3 c4 w0 - (3/5) c4^2 w0^2 - (4/15) c4 c2.

Also: sanity-check the D-formula on p = (5/4)w^4 - 3w^2 + 3/4 (constructed
to have a rescued cusp at w0=0: expect D(0) = 1 > 0, D(+-sqrt(6/5)) = -1/8).
And on F4's seed c1... (-w^4 + w^3 - 27/10 w^2 + 17/10 w): expect D < 0.
"""
import sympy as sp
from sympy import symbols, Rational as R, diff
import itertools

w = symbols("w")

def seed(c2, c3, c4):
    S1 = -1 - c2 - c3 - c4
    S2 = -(R(c2,1)/3 + R(c3,1)/4 + R(c4,1)/5)
    return [sp.nsimplify(2*S2 - S1), sp.nsimplify(2*(S1 - S2)), R(c2,1), R(c3,1), R(c4,1)]

def D_at(pco, w0):
    c0,c1,c2,c3,c4 = pco
    return sp.simplify(R(c3)**2/16 - R(3,10)*c3*c4*w0 - R(3,5)*R(c4)**2*w0**2 - R(4,15)*c4*c2)

# --- sanity seeds (not necessarily normalized) ---
for label, pco in [("constructed (5/4,0,-3,0,3/4)", [R(3,4),0,-3,0,R(5,4)]),
                   ("F4 family seed", [0, R(17,10), R(-27,10), 1, -1])]:
    p = sum(ci*w**i for i,ci in enumerate(pco))
    crits = sp.nroots(diff(p,w), n=30)
    rows = []
    for r in crits:
        if abs(sp.im(r)) < 1e-20:
            rows.append((float(sp.re(r)), float(sp.N(D_at(pco, sp.Float(sp.re(r),40)), 15))))
    print(f"{label}: cusps (w0, D) = {rows}")

# --- stats over the normalized grid ---
from collections import Counter
stats = Counter(); solo = Counter(); rescue_by_count = Counter()
for c2, c3, c4 in itertools.product(range(-4,5), range(-4,5), range(-4,5)):
    if c4 == 0: continue
    pco = seed(c2,c3,c4)
    p = sum(ci*w**i for i,ci in enumerate(pco))
    kap = diff(p, w).subs(w, 1)
    if kap == -2: continue
    crits = [r for r in sp.nroots(diff(p,w), n=30, maxsteps=200) if abs(sp.im(r)) < 1e-20]
    n = len(crits)
    Ds = [sp.N(D_at(pco, sp.Float(sp.re(r),40)), 15) for r in crits]
    signs = tuple("+" if D > 0 else ("-" if D < 0 else "0") for D in Ds)
    stats[(n, signs)] += 1
    solo[n] += 1
print("\ncusp-count distribution:", dict(solo))
print("(n_cusps, D-sign pattern) -> #seeds:")
for k in sorted(stats, reverse=True):
    print("   ", k, stats[k])
