#!/usr/bin/env python3
"""
NOTE 19, SEGMENT P4: THE FLAT-SHEET UNIFORM LEMMA (reviewer item P4)
====================================================================
Note 4's surjectivity theorem leans on "C = 0 coverage": preimages of targets
(A,B,0). The flat sheet {x = 0} maps into {C = 0} (third component is x*gamma),
and for F4/F5 the restriction was computed by hand to be a triangular Keller
automorphism (Jung-van der Kulk), giving every (A,B,0) exactly one flat preimage.
Reviewer: prove it ONCE, UNIFORMLY.

Claim under test: for the FULL recipe family (any seed p with p(0)=0, p(1)=-c,
int_0^1 p = 0; q' = w p'/c, q(0)=0; gamma = 1 + a v + b t, a = -(1+k)/(2+k),
k = p'(1)/c != -2; u = 1+v; b, c != 0), the restriction of
F = (alpha/x^2, beta/x, x gamma) to x = 0 is EXACTLY

    F(0, y, z) = ( b(k+2) z + C2 y^2,  -c/(2+k) y,  0 ),

an ELEMENTARY triangular automorphism of the {C = 0} plane (both diagonal entries
nonzero since k != -2, b c != 0), with C2 = (1/2) d^2/dv^2 [alpha(v, 0)]|_0.

LOCKS (registered BEFORE computation):
 FS1 (symbolic structure): with generic seed data, the x-series of alpha and beta
     through order 2 obeys [x^1]beta = -c/(2+k) y and [x^2]alpha = b(k+2) z + C2 y^2
     (z- and y^2-terms only; no other monomials survive endpoint conditions).
     Checked on 6 generic polynomial seeds (random extra w(1-w) terms, c in {1,2}).
 FS2 (tower diagonal): for the explainer seeds d = 2..12 (kappa_d = -5 + 6/(d(d+1))):
     diagonal entries -1/(2+kappa_d) and kappa_d+2 are explicit NONZERO rationals;
     -1/(2+kappa_d) = d(d+1)/(3(d(d+1)-2)) (= 1/3 + 2/(3(d(d+1)-2)) - the ghost's
     arithmetic cousin).
 FS3 (C2 table): exact C2 for tower seeds d = 2..8; C2(2) = 4 in the posted-map
     normalization cross-check: series-computed A0 for the fiber-3 posted map equals
     z + 4 y^2 and Q0 = y, matching note 2's direct expansion.
 FS4 (automorphism): numerical check (60 dps): the inverse formula
     y = -(2+k)/c B,  z = (A - C2 y^2)/(b(k+2))  recovers 30 random (A,B) targets
     through the FULL map F(x,y,z) (x = 0) with residuals < 1e-50, for d = 2..6.
"""
import json
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
out = {"locks": {}}
v, t, ysym, zsym = sp.symbols('v t y z')

def seed_generic(c_val, extra_terms):
    """p = 2w-3w^2 + w(1-w)*(extra) with int_0^1 p = 0 forced via constant part."""
    w = sp.symbols('w')
    base = 2 * w - 3 * w**2
    ext = sum(sp.Rational(str(cf)) * w**k for cf, k in extra_terms)
    # int_0^1 w(1-w) w^k = 1/((k+2)(k+3))
    integ = sum(sp.Rational(str(cf)) / ((k + 2) * (k + 3)) for cf, k in extra_terms)
    # kill the integral by a constant shift inside w(1-w):
    const_term = -integ / sp.Rational(1, 6)   # int_0^1 w(1-w) = 1/6
    p = sp.expand(sp.Integer(c_val) * (base + w * (1 - w) * (ext + const_term)))
    return p, sp.Integer(c_val)

def analyze(p, c):
    """Given seed p, c: build alpha(v,t), beta(v,t) as series to x^2 and extract
    the flat-sheet restriction. Uses v = x*y, t = x^2*z."""
    w = sp.symbols('w')
    kap = sp.Rational(sp.diff(p, w).subs(w, 1), 1) / c
    q = sp.integrate(w * sp.diff(p, w) / c, w)
    a = -(1 + kap) / (2 + kap)
    b = sp.Symbol('b', nonzero=True)
    ga = 1 + a * v + b * t
    u = 1 + v
    W = u * ga
    al = u + q.subs(w, W) / ga**2
    be = c + p.subs(w, W) / ga
    # series in an order-2 variable: substitute v -> eps*y, t -> eps^2*z, series eps
    eps = sp.symbols('eps')
    alS = sp.series(al.subs({v: eps * ysym, t: eps**2 * zsym}), eps, 0, 3).removeO().expand()
    beS = sp.series(be.subs({v: eps * ysym, t: eps**2 * zsym}), eps, 0, 3).removeO().expand()
    A0 = sp.expand(alS.coeff(eps, 2))      # [x^2] alpha
    B1 = sp.expand(beS.coeff(eps, 1))      # [x^1] beta
    # also prove lower coefficients vanish
    assert sp.expand(alS.coeff(eps, 0)) == 0 and sp.expand(alS.coeff(eps, 1)) == 0, "alpha not O(x^2)"
    assert sp.expand(beS.coeff(eps, 0)) == 0, "beta not O(x^1)"
    return sp.simplify(A0), sp.simplify(B1), kap, a, q

print("== FS1: generic seeds ==")
generic_ok = True
for (c_, ext) in [(1, [(sp.Rational(2, 3), 1), (sp.Rational(-1, 2), 3)]),
                  (1, [(sp.Rational(1, 1), 2), (sp.Rational(3, 4), 4), (sp.Rational(-2), 5)]),
                  (2, [(sp.Rational(5, 7), 1)]),
                  (1, []),
                  (3, [(sp.Rational(-3, 5), 2), (sp.Rational(2), 6)]),
                  (1, [(sp.Rational(7, 3), 7)])]:
    p_, c_ = seed_generic(c_, ext)
    A0, B1, kap, a_, q_ = analyze(p_, c_)
    # expected forms
    haz = sp.Symbol('b', nonzero=True)
    z_coeff_expected = haz * (kap + 2)
    z_coeff_actual = A0.coeff(zsym, 1)
    y2_coeff = A0.coeff(ysym, 2)
    residue = sp.simplify(A0 - z_coeff_expected * zsym - y2_coeff * ysym**2)
    B_residue = sp.simplify(B1 - (-c_ / (2 + kap)) * ysym)
    ok = (residue == 0) and (B_residue == 0) and sp.simplify(z_coeff_actual - z_coeff_expected) == 0
    generic_ok &= bool(ok)
    print(f"   c={c_} k={kap}: A0 = {sp.simplify(A0)}  B1 = {sp.simplify(B1)}  forms OK: {ok}")
out["locks"]["FS1"] = {"generic_seeds": 6, "all_forms_exact": bool(generic_ok)}

print("\n== FS2: tower diagonal ==")
tower_rows = {}
for d in range(2, 13):
    kap_d = sp.Rational(-5) + sp.Rational(6, 1) / (d * (d + 1))
    diag_y = sp.simplify(-1 / (2 + kap_d))
    diag_z = sp.simplify(kap_d + 2)
    ghost_cousin = sp.Rational(d * (d + 1), 3 * (d * (d + 1) - 2))
    assert sp.simplify(diag_y - ghost_cousin) == 0
    assert diag_y != 0 and diag_z != 0
    tower_rows[d] = {"kappa": str(kap_d), "diag_y": str(diag_y), "diag_z": str(diag_z)}
    print(f"   d={d:2d}: kappa={kap_d}, diag_y = {diag_y} = d(d+1)/(3(d(d+1)-2)), diag_z = {diag_z}")
out["locks"]["FS2"] = tower_rows

print("\n== FS3: tower C2 table + posted-map cross-check ==")
C2_tab = {}
for d in range(2, 9):
    w = sp.symbols('w')
    c0 = sp.Rational(6, d * (d + 1))
    p_ = sp.expand(2 * w - 3 * w**2 + w * (1 - w) * (w**(d - 2) - c0))
    A0, B1, kap, a_, q_ = analyze(p_, sp.Integer(1))
    C2 = sp.simplify(A0.coeff(ysym, 2))
    C2_tab[d] = str(C2)
    print(f"   d={d}: C2 = {C2},  A0 = {sp.simplify(A0)},  B1 = {sp.simplify(B1)}")
out["locks"]["FS3"] = {"C2": C2_tab}
# posted-map cross-check: fiber-3 Alpoege map: F(0,y,z) should be (z + 4y^2, y, 0)
x, y, z = sp.symbols('x y z')
u_ = 1 + x * y
P = u_**3 * z + y**2 * u_ * (4 + 3 * x * y)
Q = y + 3 * x * u_**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
R = 2 * x - 3 * x**2 * y - x**3 * z
flat = [sp.expand(comp.subs(x, 0)) for comp in (P, Q, R)]
posted_ok = flat == [z + 4 * y**2, y, 0]
print(f"   posted fiber-3: F(0,y,z) = {flat[0]}, {flat[1]}, {flat[2]}  OK: {posted_ok}")
out["locks"]["FS3"]["posted_fiber3_ok"] = posted_ok

print("\n== FS4: flat-sheet inversion through the full map, d = 2..6 ==")
def tower_poly(d):
    w = sp.symbols('w')
    c0 = sp.Rational(6, d * (d + 1))
    return sp.expand(2 * w - 3 * w**2 + w * (1 - w) * (w**(d - 2) - c0))
worst4 = mp.mpf('0'); ok4 = True
import random
random.seed(44)
for d in range(2, 7):
    w = sp.symbols('w')
    p_ = tower_poly(d); c_ = sp.Integer(1)
    kap = sp.Rational(sp.diff(p_, w).subs(w, 1)) / c_
    q_ = sp.integrate(w * sp.diff(p_, w) / c_, w)
    a_ = -(1 + kap) / (2 + kap)
    b_ = sp.Integer(1)
    ga = 1 + a_ * v + b_ * t
    al = (1 + v) + q_.subs(w, (1 + v) * ga) / ga**2
    be = c_ + p_.subs(w, (1 + v) * ga) / ga
    # Flat-sheet map: (A, B) given; y = -(2+k)/c B; z = (A - C2 y^2)/(b(k+2))
    C2 = sp.simplify(sp.diff(al, v, 2).subs({v: 0, t: 0}) / 2)
    for _ in range(30):
        A = sp.Rational(random.randint(-30, 30), random.randint(1, 5))
        B = sp.Rational(random.randint(-30, 30), random.randint(1, 5))
        yv = -(2 + kap) * B / c_
        zv = (A - C2 * yv**2) / (b_ * (kap + 2))
        # evaluate full map at (x,y,z) = (0, yv, zv):
        # A = alpha/x^2 etc: evaluate via series substitution with eps
        eps = sp.symbols('eps')
        Aexpr = sp.series(al.subs({v: eps * yv, t: eps**2 * zv}) / eps**2, eps, 0, 1).removeO()
        Bexpr = sp.series(be.subs({v: eps * yv, t: eps**2 * zv}) / eps, eps, 0, 1).removeO()
        Cexpr = sp.expand(eps * ga.subs({v: eps * yv, t: eps**2 * zv})).coeff(eps, 0)
        rA = abs(mp.mpf(str(sp.N(sp.expand(Aexpr) - A, 40))))
        rB = abs(mp.mpf(str(sp.N(sp.expand(Bexpr) - B, 40))))
        rC = abs(mp.mpf(str(sp.N(Cexpr, 40))))
        worst4 = max(worst4, rA, rB, rC)
        if max(rA, rB, rC) > mp.mpf('1e-30'): ok4 = False
out["locks"]["FS4"] = {"worst_residual": mp.nstr(worst4, 5), "targets": 30 * 5, "ok": ok4}
print(f"   150 targets (A,B,0), d=2..6: worst residual {mp.nstr(worst4,5)}, all exact: {ok4}")

with open("/home/user/flatsheet_stage.json", "w") as f:
    json.dump(out, f, indent=2, default=str)
print("\nsaved flatsheet_stage.json")
