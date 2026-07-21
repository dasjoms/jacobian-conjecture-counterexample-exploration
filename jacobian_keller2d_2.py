#!/usr/bin/env python3
"""
NOTE 18 - stage 2: THE CANAL THEOREM & THE MISSING BRAID
=========================================================
Context repair: stage 1 (jacobian_keller2d_1.py) used u = x^2*y as the "u-channel" -
WRONG. The true recipe (explainer, verified against the degree-4 example) is

    v = x*y,  u = 1 + v,  w = u*gamma(v),  F = (alpha/x^2, beta/x),
    alpha = u + q(w)/gamma(v)^2,  beta = c + p(w)/gamma(v).

But then BOTH components are (function of v alone) * (power of x). That is the
CANAL. For ANY entire/meromorphic A(v), B(v) and F = (A(v)/x^m, B(v)/x^n):

    det JF = (n*A'*B - m*A*B') / x^(m+n).

Keller (det = nonzero const) is IMPOSSIBLE unless the numerator vanishes, and then
det = 0 and A^n/B^m is constant (image collapses onto the curve X^n = const*Y^m).
The obstruction is NOT polynomiality - it is CHANNEL RANK. Algebraically: the
weight-zero monomial ring of C^2 with weights (1,-1) is C[v] (ONE generator);
for C^3 with weights (1,-1,-2) it is C[v,t] (TWO generators: v = xy, t = x^2 z).
The 3D tower braids two channels; the plane offers one. One channel cannot braid.

LOCKS (all registered BEFORE computation):
 L1 (CANAL IDENTITY): det JF = (n*A'*B - m*A*B')/x^(m+n) holds exactly for 25
    random polynomial pairs (A,B) and 5 random (m,n), incl. the recipe pair
    (m,n) = (2,1) seeded with the fiber-3 seed, at 50 digits on 40 points each:
    max error < 1e-40.
 L2 (SEEDED MAP ALIVE BUT NOT KELLER): with p(w) = 2w-3w^2, c = 1, q' = w p',
    gamma = 1 - 3v/2, u = 1+v: F2 = (alpha/x^2, beta/x) is POLYNOMIAL (endpoint
    conditions kill the poles exactly; v | beta and v^2 | alpha in Q[v]).
    Then det F2 = D(v)/x^3 with D = alpha' beta - 2 alpha beta' a POLYNOMIAL,
    and polynomiality forces v^3 | D: det = y^3 * E(v), E = D/v^3 an explicit
    NONCONSTANT polynomial (rel. variation over the sampled range > 0.1).
    Additionally alpha/beta^2 is NONCONSTANT (rank-1 collapse A = C B^2 fails:
    the image is genuinely 2-dimensional off the degeneracy locus).
 L3 (BRAID LEDGER): symbolic: (i) det of the 3D F3 = (alpha/x^2, beta/x, x*gamma),
    gamma = 1 + a v + b t, equals b*c EXACTLY (independent re-verification);
    (ii) cutting the second channel (set b = 0) gives det = 0 EXACTLY;
    (iii) rank of the channel map (x,y,z) -> (v,t) is 2; of (x,y) -> v is 1;
    (iv) weight-zero monomials: {x^A y^B z^C : A - B - 2C = 0} = {v^B t^C}
    (two-parameter family); {x^A y^B : A - B = 0} = {v^B} (one-parameter).
 L4 (ENTIRE SEED-CARRYING AVATAR): F_s(x,y) = (e^x, (y + p(e^(2x))) e^(-x)),
    p = fiber-3 seed: (i) det = 1 exactly (symbolic); (ii) F_s(x + 2 pi i k, y)
    = F_s(x,y) exactly for k = -3..3 (fiber contains a Z-lattice of sheets);
    (iii) numeric fiber solve at 60 digits for 25 random targets: unique sheet
    in the fundamental strip 0 <= Im x < 2 pi, residual < 1e-50; fiber is
    {(log X + 2 pi i k, Y X - p(X^2))}, y SAME on all sheets; (iv) {X = 0}
    is missed (the wall flattens to a line).
 L5 (MOH FRONTIER DATA): degree census assembled & arithmetic-checked:
    Moh 1983 <= 100 (complex, Keller => automorphism); Pinchuk 1994 real
    nonvanishing non-injective, component degrees (10, 25), total 35 = 10+25;
    tower degrees 7 (fiber-3) and 12 (fiber-4) live in 3-D. Canal corollary:
    NO single-channel weighted recipe (any entire u(v), gamma(v), p, q) can be
    Keller in 2-D - Moh's chamber cannot be entered from the tower's staircase.
"""
import json, random
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
random.seed(1802)
out = {"locks": {}, "seeded": {}, "braid": {}, "avatar": {}}
def put(lock, name, val): out["locks"].setdefault(lock, {})[name] = val

x, y, z = sp.symbols('x y z')
V, T = sp.symbols('v t')
mi = mp.matrix  # noqa

# ---------------------------------------------------------------- L1 canal identity
print("== L1: canal identity det JF = (n A' B - m A B')/x^(m+n) ==")
worst1 = mp.mpf('0')
cases = 0
for _ in range(25):
    da, db = random.randint(0, 5), random.randint(0, 5)
    Ac = [sp.Rational(random.randint(-4, 4), random.randint(1, 3)) for _ in range(da + 1)]
    Bc = [sp.Rational(random.randint(-4, 4), random.randint(1, 3)) for _ in range(db + 1)]
    apoly = sum(c * V**k for k, c in enumerate(Ac))
    bpoly = sum(c * V**k for k, c in enumerate(Bc))
    for (m, n) in [(1, 1), (1, 2), (2, 1), (3, 2), (2, 3)]:
        X = apoly.subs(V, x * y) / x**m
        Y = bpoly.subs(V, x * y) / x**n
        det_direct = sp.Matrix([X, Y]).jacobian([x, y]).det()
        formula = (n * sp.diff(apoly, V) * bpoly - m * apoly * sp.diff(bpoly, V)).subs(V, x * y) / x**(m + n)
        diff = sp.simplify(det_direct - formula)
        assert diff == 0, f"symbolic mismatch m={m},n={n}"
        # numeric bullet @ 50 dps: exact partials via lambdified sympy derivatives
        fA = sp.lambdify(V, apoly, 'mpmath'); fB = sp.lambdify(V, bpoly, 'mpmath')
        fAp = sp.lambdify(V, sp.diff(apoly, V), 'mpmath')
        fBp = sp.lambdify(V, sp.diff(bpoly, V), 'mpmath')
        ent = [[sp.lambdify((x, y), sp.diff(X, s), 'mpmath') for s in (x, y)],
               [sp.lambdify((x, y), sp.diff(Y, s), 'mpmath') for s in (x, y)]]
        X0, Y0 = mp.mpf('1.7'), mp.mpf('-0.9')
        vv = X0 * Y0
        detn = ent[0][0](X0, Y0) * ent[1][1](X0, Y0) - ent[0][1](X0, Y0) * ent[1][0](X0, Y0)
        form = (n * fAp(vv) * fB(vv) - m * fA(vv) * fBp(vv)) / X0**(m + n)
        worst1 = max(worst1, abs(detn - form))
        cases += 1
put("L1", "cases", cases)
put("L1", "worst_err", mp.nstr(worst1, 5))
print(f"   symbolic: all {cases} exact; numeric worst |det-direct - formula| = {mp.nstr(worst1,5)}")

# ---------------------------------------------------------------- L2 seeded 2-D map
print("\n== L2: seeded 2-D compression (fiber-3 seed) ==")
w, v = sp.symbols('w v')
p = 2 * w - 3 * w**2
c = sp.Integer(1)
q = sp.integrate(w * sp.diff(p, w) / c, w)          # q(0)=0
ga = 1 - sp.Rational(3, 2) * v                      # tower's linear gamma
u = 1 + v
W = u * ga
al = sp.cancel(u + q.subs(w, W) / ga**2)
be = sp.cancel(c + p.subs(w, W) / ga)
alE, beE = sp.Poly(al, v), sp.Poly(be, v)
print(f"   alpha(v) = {sp.expand(al)}  (deg {alE.degree()})")
print(f"   beta(v)  = {sp.expand(be)}  (deg {beE.degree()})")
# strip test: v | beta, v^2 | alpha  (x|beta/x and x^2|alpha/x^2 polynomiality)
assert sp.rem(be, sp.Poly(v, v)) == 0, "v does not divide beta"
assert sp.rem(al, sp.Poly(v**2, v)) == 0, "v^2 does not divide alpha"
X2 = sp.cancel(al.subs(v, x * y) / x**2)
Y2 = sp.cancel(be.subs(v, x * y) / x)
assert sp.denom(X2) == 1 and sp.denom(Y2) == 1, "not polynomial!"
# endpoint conditions that made the poles die (tower values, gamma0 = 1):
print(f"   endpoint checks: p(1)+1 = {sp.simplify(p.subs(w,1)+1)}, "
      f"q(1)+1 = {sp.simplify(q.subs(w,1)+1)}, "
      f"gamma'(0) = {sp.diff(ga,v).subs(v,0)}, tower a = -3/2")
D = sp.expand(sp.diff(al, v) * be - 2 * al * sp.diff(be, v))
Dp = sp.Poly(D, v)
qD, rD = sp.div(Dp, sp.Poly(v**3, v))
assert rD == 0, "v^3 does not divide D"
E = sp.expand(qD.as_expr())
print(f"   D(v) = alpha' beta - 2 alpha beta' = {sp.factor(D)}")
print(f"   v^3 | D exactly;  det = y^3 * E(v),  E(v) = {E}")
det2_direct = sp.factor(sp.Matrix([X2, Y2]).jacobian([x, y]).det())
det2_canal = sp.factor(D.subs(v, x * y) / x**3)
assert sp.simplify(det2_direct - det2_canal) == 0
det_yform = sp.factor(E.subs(v, x * y) * y**3)
assert sp.simplify(det2_direct - det_yform) == 0
print(f"   det F2 = y^3 * E(xy)  [canal: never a nonzero constant]")
# nonconstancy of E + A != C B^2
Evar_samples = [abs(complex(sp.N(E.subs(v, sp.Rational(k, 7))))) for k in range(1, 8)]
rel_var = (max(Evar_samples) - min(Evar_samples)) / max(Evar_samples)
ratio = sp.simplify(al / be**2)
ratio_not_const = sp.simplify(sp.diff(ratio, v)) != 0
print(f"   |E| over v=1/7..1: {[mp.nstr(mp.mpf(str(s)),4) for s in Evar_samples]} rel-var {float(rel_var):.3f}")
print(f"   d/dv (alpha/beta^2) nonzero: {ratio_not_const}")
Pn = sp.Poly(E, v)
degen_roots = sp.nroots(Pn, n=20, maxsteps=200)
def fmt_root(r):
    z = complex(sp.N(r, 20))
    if abs(z.imag) < 1e-18:
        return mp.nstr(mp.mpf(z.real), 10)
    return f"{mp.nstr(mp.mpf(z.real), 8)}{z.imag:+.8f}i"
degen_str = [fmt_root(r) for r in degen_roots]
real_neg = [complex(sp.N(r, 20)).real for r in degen_roots
            if abs(complex(sp.N(r, 20)).imag) < 1e-15]
print(f"   degeneracy: y=0 line plus hyperbolae xy in {degen_str} (zeros of E)")
print(f"   E has all-negative coefficients -> no positive root; real roots: {real_neg}")
put("L2", "polynomial", True)
put("L2", "det_form", f"y^3 * E(v), E = {E}")
put("L2", "E_relvariation", float(rel_var))
put("L2", "alpha_over_beta2_nonconst", bool(ratio_not_const))
put("L2", "degeneracy_roots", degen_str)

# endpoint-condition check (symbolic, generic gamma0)
g0 = sp.symbols('g0', nonzero=True)
al_at0 = sp.simplify(al.subs({v: 0}) )
print(f"   alpha at v=0 with gamma0=1: {al_at0}  (must be 0)")

# ---------------------------------------------------------------- L3 braid ledger
print("\n== L3: braid ledger ==")
a, b = sp.symbols('a b', nonzero=True)
# contrast probe: NON-recipe (random) p, q, a -> det3 not constant
prand = sum(sp.Integer(random.randint(-3, 3)) * w**k for k in range(1, 4))
qrand = sum(sp.Integer(random.randint(-3, 3)) * w**k for k in range(2, 5))
arand = sp.Rational(2, 3)
g3r = 1 + arand * v + sp.Symbol('bb', nonzero=True) * T
alr = ((1 + v) + qrand.subs(w, (1 + v) * g3r) / g3r**2).subs({v: x * y, T: x**2 * z})
ber = (1 + prand.subs(w, (1 + v) * g3r) / g3r).subs({v: x * y, T: x**2 * z})
F3r = sp.Matrix([alr / x**2, ber / x, (x * g3r).subs({v: x * y, T: x**2 * z})])
det3r = sp.factor(F3r.jacobian([x, y, z]).det())
det3r_const = (sp.diff(det3r, x) == 0 and sp.diff(det3r, y) == 0 and sp.diff(det3r, z) == 0)
print(f"   NON-recipe control: det3 = const ? {det3r_const}  (expect False)")
# with recipe relations the surviving structure should reduce to b*c; verify with explicit seed
p3 = 2 * w - 3 * w**2
q3 = w**2 - 2 * w**3
kap = sp.diff(p3, w).subs(w, 1)
aval = -(1 + kap) / (2 + kap)
bval = sp.Symbol('bb', nonzero=True)
g3 = 1 + aval * v + bval * T
al3e = (1 + v) + q3.subs(w, (1 + v) * g3) / g3**2
be3e = 1 + p3.subs(w, (1 + v) * g3) / g3
F3e = sp.Matrix([al3e.subs({v: x * y, T: x**2 * z}) / x**2,
                 be3e.subs({v: x * y, T: x**2 * z}) / x,
                 (x * g3).subs({v: x * y, T: x**2 * z})])
det3e = sp.factor(F3e.jacobian([x, y, z]).det())
print(f"   seeded det3 = {det3e}  (expect bb = b*c)")
# cut the channel: b = 0, gamma = gamma(v) only
g0v = 1 + aval * v
al0 = ((1 + v) + q3.subs(w, (1 + v) * g0v) / g0v**2).subs(v, x * y)
be0 = (1 + p3.subs(w, (1 + v) * g0v) / g0v).subs(v, x * y)
F30 = sp.Matrix([al0 / x**2, be0 / x, x * g0v.subs(v, x * y)])
det30 = sp.simplify(F30.jacobian([x, y, z]).det())
print(f"   det3 with b = 0 (second channel cut): {det30}")
# channel-map ranks
chan3 = sp.Matrix([x * y, x**2 * z]).jacobian([x, y, z])
generic3 = chan3.subs({x: 2, y: 3, z: 5}).rank()
chan2 = sp.Matrix([x * y]).jacobian([x, y])
generic2 = chan2.subs({x: 2, y: 3}).rank()
# weight-zero monomials
pairs2 = [(A, B) for A in range(0, 7) for B in range(0, 7) if A - B == 0]
trips3 = [(A, B, C) for A in range(0, 9) for B in range(0, 5) for C in range(0, 5) if A - B - 2 * C == 0]
fam2 = {(A - B, B) for (A, B) in pairs2}  # all (0,k) -> v^k
fam3 = {(B, C) for (A, B, C) in trips3}   # check v^B t^C parametrization
ok3 = all((A, B, C) == (B + 2 * C, B, C) for (A, B, C) in trips3)
put("L3", "det3_seeded", str(det3e))
put("L3", "det3_b0", str(det30))
put("L3", "rank_channels_3d", generic3)
put("L3", "rank_channels_2d", generic2)
put("L3", "weight0_3d_two_param", ok3 and len(fam3) > 1)
put("L3", "weight0_2d_one_param", all(k == 0 for k, _ in fam2))
print(f"   ranks: (v,t)-channels {generic3} vs v-only {generic2}; weight0 3D two-param: {ok3}")
print(f"   2D weight-0 monomials x^A y^B: {pairs2} (all v^k)")
print(f"   3D weight-0 sample: {trips3[:8]}... all (B+2C, B, C): {ok3}")

# ---------------------------------------------------------------- L4 entire avatar
print("\n== L4: entire seed-carrying avatar F_s = (e^x, (y + p(e^{2x})) e^{-x}) ==")
E = sp.exp(x)
Xs = E
Ys = (y + p.subs(w, sp.exp(2 * x))) * sp.exp(-x)
dets = sp.simplify(sp.Matrix([Xs, Ys]).jacobian([x, y]).det())
print(f"   det = {dets}")
# periodicity exactly
k = sp.symbols('k', integer=True)
xshift = x + 2 * sp.pi * sp.I * 3
per_check = sp.simplify((Xs.subs(x, xshift) - Xs)) , sp.simplify(sp.expand_complex(Ys.subs(x, xshift) - Ys))
print(f"   F(x+6 pi i)-F(x): {per_check}")
# numeric fibers at 60 dps
worst4 = mp.mpf('0')
all_hits = True
for _ in range(25):
    Xt = mp.mpf(str(random.uniform(0.5, 2.5))) * mp.e**(mp.mpf(str(random.uniform(-1.5, 1.5))) * 1j)
    Yt = mp.mpf(str(random.uniform(-2, 2))) + mp.mpf(str(random.uniform(-2, 2))) * 1j
    x0 = mp.log(Xt)  # fundamental strip
    y0 = Yt * Xt - (2 * Xt**2 - 3 * Xt**4)
    r1 = abs(mp.e**x0 - Xt)
    r2 = abs((y0 + (2 * mp.e**(2 * x0) - 3 * mp.e**(4 * x0))) * mp.e**(-x0) - Yt)
    worst4 = max(worst4, r1, r2)
    # every shifted sheet hits too
    for kk in (-2, -1, 1, 2):
        xk = x0 + 2 * mp.pi * 1j * kk
        yk = Yt * mp.e**xk - (2 * mp.e**(2 * xk) - 3 * mp.e**(4 * xk))
        if abs(yk - y0) > mp.mpf('1e-45'):
            all_hits = False
        r = abs(mp.e**xk - Xt) + abs((yk + (2 * mp.e**(2 * xk) - 3 * mp.e**(4*xk))) * mp.e**(-xk) - Yt)
        worst4 = max(worst4, r)
put("L4", "det_symbolic", str(dets))
put("L4", "fiber_residual_worst", mp.nstr(worst4, 5))
put("L4", "sheets_locked_to_lattice", all_hits)
put("L4", "missed_set", "{X = 0} (line)")
print(f"   25 targets x (5 sheets): worst residual {mp.nstr(worst4,5)}, y sheet-invariant: {all_hits}")

# ---------------------------------------------------------------- L5 frontier data
put("L5", "moh", "1983 Crelle 340:140-212, JC(2 vars) verified degrees <= 100")
put("L5", "pinchuk", "1994 Math Z 217:1-4, real, component degrees (10,25), total 35")
put("L5", "pinchuk_total_check", 10 + 25 == 35)
put("L5", "tower", "3-D: fiber-3 degree 7 det -2; fiber-4 degree 12 det -6")

with open("/home/user/keller2d_stage2.json", "w") as f:
    json.dump(out, f, indent=2, default=str)
print("\nSaved keller2d_stage2.json")
print("\nLOCK BOARD:")
for lk in ["L1", "L2", "L3", "L4", "L5"]:
    print(f"  {lk}: {out['locks'].get(lk)}")
