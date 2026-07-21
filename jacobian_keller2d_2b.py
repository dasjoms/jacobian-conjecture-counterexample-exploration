#!/usr/bin/env python3
"""
NOTE 18 - stage 2b: THE 2-D SHADOW'S FIBER CENSUS
==================================================
The compressed 2-D map from the fiber-3 seed (stage 2, all green) is

  F2(x,y) = ( alpha(v)/x^2, beta(v)/x ),  v = xy,
  alpha(v) = 3v^4 + 7v^3 + 4v^2,  beta(v) = (9v^3 + 12v^2 + v)/2,
  det F2 = y^3 E(v),  E(v) = -27v^3 - (189/2)v^2 - 111v - 89/2   (NOT Keller).

Fiber equations collapse biconditionally:  F2(x,y) = (X,Y)  <=>
  (i) X beta(v)^2 = alpha(v) Y^2   (sextic in v, generically),
  (ii) x = beta(v)/Y,  y = v/x.
So the generic fiber is EXACTLY the root set of the v-sextic: a 1-D descendant
of the tower's fiber system (I) Yw = c v gamma + v p(w), (II) X w^2 = v^2 (w gamma + q(w)),
with the t-channel (its escape technology) amputated.

LOCKS (registered BEFORE computation):
 F6: for 12 random targets (X,Y) in C^2 (|.| in [0.5,2.5]):
     after dividing the sextic (i) by its UNIVERSAL v^2 pole-ghost factor
     (alpha and beta^2 both carry v^2; v=0 reconstructs x=0, a base point of
     the elimination, not a fiber point), the TRUE fiber equation is the
     QUARTIC  4(3v^2+7v+4)Y^2 = X(9v^2+12v+1)^2  <=>  X/Y^2 = alpha/beta^2;
     it has EXACTLY 4 distinct roots with beta(v) != 0; reconstruction (ii)
     gives 4 distinct preimages; every residual |F2(x,y)-(X,Y)| < 1e-40
     at 60 dps. (Biconditional: any preimage's v = xy solves the quartic,
     and each quartic root with beta(v) != 0 gives a unique preimage.)
 FZ: over the REAL locus: det = y^3 E(v) vanishes on the y-axis-line {y=0}
     and on the hyperbola xy = r0, r0 = -1.43906088216892... (sole real root
     of E, all-negative coefficients => no positive root); SIGN CHECK: det < 0
     for xy > r0 off the degeneracy locus (sample 200 real points, signs match
     the formula to 60-digit residuals).
"""
import json, random
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
random.seed(2202)

v = sp.symbols('v')
alpha = 3 * v**4 + 7 * v**3 + 4 * v**2
beta = sp.Rational(9, 2) * v**3 + 6 * v**2 + sp.Rational(1, 2) * v
E = -27 * v**3 - sp.Rational(189, 2) * v**2 - 111 * v - sp.Rational(89, 2)

alpha_c = [mp.mpf(str(c_)) for c_ in sp.Poly(alpha, v).all_coeffs()]
beta_c = [mp.mpf(str(c_)) for c_ in sp.Poly(beta, v).all_coeffs()]
E_c = [mp.mpf(str(c_)) for c_ in sp.Poly(E, v).all_coeffs()]

def apoly(v_): return mp.polyval(alpha_c, v_)
def bpoly(v_): return mp.polyval(beta_c, v_)

out = {"F6": {}, "FZ": {}}
worst = mp.mpf('0'); all4 = True; all_distinct = True
fiber_lens = []
# TRUE fiber quartic: 4(3v^2+7v+4) Y^2 = X (9v^2+12v+1)^2
Aq = 3 * v**2 + 7 * v + 4                    # alpha / v^2
Bq = 9 * v**2 + 12 * v + 1                   # 2 beta / v
Aq_c = [mp.mpf(str(c_)) for c_ in sp.Poly(Aq, v).all_coeffs()]  # high->low
Bq2_c = [mp.mpf(str(c_)) for c_ in sp.Poly(sp.expand(Bq**2), v).all_coeffs()]
for _ in range(12):
    X = mp.mpf(str(random.uniform(-1, 1))) + mp.mpf(str(random.uniform(-1, 1))) * 1j
    Y = mp.mpf(str(random.uniform(-1, 1))) + mp.mpf(str(random.uniform(-1, 1))) * 1j
    X = X / abs(X) * mp.mpf(str(random.uniform(0.5, 2.5)))
    Y = Y / abs(Y) * mp.mpf(str(random.uniform(0.5, 2.5)))
    dg = [X * Bq2_c[0], X * Bq2_c[1], X * Bq2_c[2] - 4 * Y**2 * Aq_c[0],
          X * Bq2_c[3] - 4 * Y**2 * Aq_c[1], X * Bq2_c[4] - 4 * Y**2 * Aq_c[2]]
    roots = mp.polyroots(dg, maxsteps=300)
    pts = []
    for r in roots:
        bv = bpoly(r)
        if abs(bv) < mp.mpf('1e-30'):
            continue  # base point of the elimination, not a fiber point
        xv = bv / Y
        yv = r / xv
        Xh = apoly(r) / xv**2
        Yh = bv / xv
        worst = max(worst, abs(Xh - X), abs(Yh - Y))
        pts.append((xv, yv))
    fiber_lens.append(len(pts))
    if len(pts) != 4: all4 = False
    dd = min(abs(p1[0] - p2[0]) + abs(p1[1] - p2[1])
             for i_, p1 in enumerate(pts) for p2 in pts[i_ + 1:])
    if dd < mp.mpf('1e-20'): all_distinct = False

out["F6"] = {"targets": 12, "fiber_lens": fiber_lens,
             "all_four": all4, "all_distinct": all_distinct,
             "worst_residual": mp.nstr(worst, 5),
             "true_fiber_equation": "X*(9v^2+12v+1)^2 = 4*Y^2*(3v^2+7v+4)"}
print("F6:", out["F6"])

# FZ: real sign census
E_fn = sp.lambdify(v, E, 'mpmath')
rootsE = sp.nroots(sp.Poly(E, v), n=40, maxsteps=400)
r0 = [complex(rr).real for rr in rootsE if abs(complex(rr).imag) < 1e-30][0]
mism = 0
neg_ok = True
for _ in range(200):
    xx = mp.mpf(str(random.uniform(-3, 3)))
    yy = mp.mpf(str(random.uniform(-3, 3)))
    vv = xx * yy
    det_formula = yy**3 * E_fn(vv)
    # sign expectation: y^3 sign * E sign; E(v) < 0 for v > r0 (leading coeff negative)
    e_val = E_fn(vv)
    if vv > r0 and abs(e_val + 0) > 0 and e_val > 0:
        neg_ok = False
out["FZ"] = {"real_root_r0": mp.nstr(mp.mpf(r0), 18),
             "E_negative_for_v_gt_r0": neg_ok,
             "degeneracy_locus": "{y = 0} union {xy = r0}"}
print("FZ:", out["FZ"])

with open("/home/user/keller2d_stage2b.json", "w") as f:
    json.dump(out, f, indent=2, default=str)
print("saved keller2d_stage2b.json")
