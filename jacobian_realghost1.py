"""
Note 5, stage 1: the d=4 surjective Keller map F4 — exact fiber machinery,
C=0 frontier, real-cusp decider, and rational certificates.

Conventions (c = b = 1):
  gamma = 1 + a*v + t,  v = x*y,  t = x^2*z,  u = 1 + v,  w = u*gamma
  F = (alpha/x^2, beta/x, x*gamma), alpha = u + q(w)/gamma^2, beta = 1 + p(w)/gamma
Fiber lemma (C != 0): preimages <-> roots w* of  h(w) := Phi(w) - BC*w + AC^2 = 0
  with gamma(w*) := BC - p(w*) != 0, via
    x = C/g, u = w*/g, y = (u-1)/x, t = g - 1 - a*(u-1), z = t/x^2 .
Phi(w) = int_0^w p.
"""
import sympy as sp
from sympy import Rational as R, symbols, expand, factor, cancel, simplify, Matrix, Poly, nroots, roots

x, y, z, w = symbols("x y z w")

def build(d):
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))
    q = sp.expand(sp.integrate(w*sp.diff(p, w), w))
    kap = sp.diff(p, w).subs(w, 1)
    a = sp.Rational(-(1 + kap) / (2 + kap))
    u = 1 + x*y
    g = 1 + a*x*y + x**2*z
    ws = u*g
    alpha = u + q.subs(w, ws)/g**2
    beta = 1 + p.subs(w, ws)/g
    f1 = sp.expand(sp.cancel(alpha/x**2))
    f2 = sp.expand(sp.cancel(beta/x))
    f3 = sp.expand(x*g)
    return sp.expand(p), sp.expand(sp.integrate(p, w)), q, a, u, g, f1, f2, f3

p, Phi, q, a_coef, u_ex, g_ex, f1, f2, f3 = build(4)
assert f1.is_polynomial(x, y, z) and f2.is_polynomial(x, y, z)
print("p4(w) =", p)
print("Phi4(w) =", Phi, "= int_0^w p4")
print("q4(w) =", q, "  a =", a_coef, "  kappa =", sp.diff(p, w).subs(w, 1))
assert sp.expand(Phi - (w*p - q)) == 0, "Phi = wp - q"
print("identity Phi = w*p - q: OK")

# ---------------- C = 0 frontier ----------------
print("\n=== C = 0: flat sheet x = 0 ===")
S1 = sp.expand(f1.subs(x, 0)); S2 = sp.expand(f2.subs(x, 0))
print("S(y,z) = (f1(0,y,z), f2(0,y,z)) =")
print("  S1 =", S1)
print("  S2 =", S2)
JS = Matrix([S1, S2]).jacobian([y, z]).det()
print("  det J(S) =", sp.factor(JS))

print("\n=== C = 0: gamma-sheet (x != 0, gamma = 0) ===")
# parametrize gamma=0 surface: z = -(1 + a*x*y)/x^2 ; restrict f1, f2
zs = -(1 + a_coef*x*y)/x**2
f1g = sp.simplify(sp.expand(f1.subs(z, zs))*x**30)  # clear denom crudely then re-divide
f2g = sp.simplify(sp.expand(f2.subs(z, zs))*x**30)
# better: substitute and cancel properly
f1g = sp.cancel(f1.subs(z, zs))
f2g = sp.cancel(f2.subs(z, zs))
print("  f1|_{gamma=0} =", sp.expand(f1g))
print("  f2|_{gamma=0} =", sp.expand(f2g))
uu = 1 + x*y
p1 = sp.Rational(17, 10); q2 = sp.Rational(17, 20)   # p'(0), q'' (2!*coeff of w^2)
expect1 = sp.expand(uu*(1 + q2*uu)/x**2)
expect2 = sp.expand((1 + p1*uu)/x)
print("  matches  f1 = u(1+q2*u)/x^2   :", sp.simplify(f1g - expect1) == 0)
print("  matches  f2 = (1 + p1*u)/x    :", sp.simplify(f2g - expect2) == 0)

# gamma-sheet fiber over (A,B,0): x = (1+p1*u)/B  and  A*x^2 = u(1+q2*u)
#   => (1+p1*u)^2 * A = B^2 * u * (1+q2*u)  : quadratic in u => generically 2 points
print("\n  gamma-sheet fiber-count over generic (A,B,0): 2 points (quadratic in u)")
# numeric spot check: target (2,3,0): solve the quadratic, verify both points map there
import math
A0, B0 = 2.0, 3.0
# A*(1+p1*u)^2 - B^2*u*(1+q2*u) = 0
c2 = A0*float(p1)**2 - B0**2*float(q2)
c1 = 2*A0*float(p1) - B0**2
c0 = A0
duu = math.sqrt(c1*c1 - 4*c2*c0)
f1n = sp.lambdify((x,y,z), f1, "math"); f2n = sp.lambdify((x,y,z), f2, "math"); f3n = sp.lambdify((x,y,z), f3, "math")
for sgn in (1, -1):
    uv = (-c1 + sgn*duu)/(2*c2)
    xv = (1 + float(p1)*uv)/B0
    yv = (uv-1)/xv
    zv = -(1+float(a_coef)*xv*yv)/xv**2
    print(f"   u={uv:.6f} x={xv:.6f} y={yv:.6f} z={zv:.6f} ->",
          f"({f1n(xv,yv,zv):.10f},{f2n(xv,yv,zv):.10f},{f3n(xv,yv,zv):.10f})")
# over (0,B,0): two gamma points x = 1/B (u=0) and x = -p1/B?? recheck: u=-1/q2: f2=(1+p1*(-1/q2))/x
print("  over (0,B,0): u=0 -> x=1/B ; u=-20/17 -> f2=(1- (17/10)(20/17))/x=-1/x=B -> x=-1/B  ... verified below")

# infinite fibers: points (0,1,0) and (0,-1,0)?
print("\n  claimed 1-dim fiber over (0,1,0): xy = -1, z = (a-1)/x^2")
pf1 = sp.simplify(f1.subs({y: -1/x, z: (a_coef-1)/x**2}))
pf2 = sp.simplify(f2.subs({y: -1/x, z: (a_coef-1)/x**2}))
pf3 = sp.simplify(f3.subs({y: -1/x, z: (a_coef-1)/x**2}))
print("   image curve =", pf1, pf2, pf3)
u0 = -1/sp.Rational(q2) if not isinstance(q2, sp.Rational) else -1/q2
print(f"  claimed 1-dim fiber over (0,-1,0): u = -1/q2 = {u0}")
yv0 = sp.simplify(u0 - 1)  # at x = 1
zv0 = sp.simplify(-(1 + a_coef*yv0))
im1 = sp.simplify(f1.subs({x: 1, y: yv0, z: zv0}))
im2 = sp.simplify(f2.subs({x: 1, y: yv0, z: zv0}))
im3 = sp.simplify(f3.subs({x: 1, y: yv0, z: zv0}))
print(f"   sample point (1, {yv0}, {zv0}) -> ({im1}, {im2}, {im3})")

# ---------------- the real-cusp decider ----------------
print("\n=== REAL cusp check (the surjectivity-over-R decider) ===")
pp = sp.diff(p, w)
print("p4'(w) =", sp.expand(pp))
real_crits = [r for r in sp.nroots(pp, n=30, maxsteps=100) if abs(sp.im(r)) < sp.Float('1e-25')]
print("real critical points of p (order-3 cusp locations):", [float(r) for r in real_crits])
for r in real_crits:
    s0 = p.subs(w, r)                # = BC   (take C = 1 -> B = s0)
    r0 = r*s0 - Phi.subs(w, r)       # = AC^2 (take C = 1 -> A = r0)
    h = sp.Poly(sp.expand(Phi - s0*w + r0).evalf(25), w)
    allr = sp.nroots(h, n=25, maxsteps=200)
    trip = [rr for rr in allr if abs(rr - r) < sp.Float('1e-12')]
    rest = [rr for rr in allr if abs(rr - r) >= sp.Float('1e-12')]
    print(f"  cusp w0 = {float(r):.10f}")
    print(f"    cusp target (A,B,C) = ({float(r0):.8f}, {float(s0):.8f}, 1)  [real target]")
    print(f"    roots of h: {[complex(rr) for rr in allr]}")
    print(f"    root cluster at w0: {len(trip)} (expect 3) ; other roots: {[complex(rr) for rr in rest]}")
    nreal = sum(1 for rr in rest if abs(sp.im(rr)) < 1e-12)
    print(f"    REAL simple partners: {nreal}  -> real preimage exists at this cusp: {nreal > 0}")
# also: are there cusps of order 4 (p' = p'' = 0 simultaneously)?
res_pp = sp.resultant(sp.diff(p, w), sp.diff(p, w, 2), w)
print("  resultant(p', p'') =", res_pp, " (nonzero => no order-4 cusp)")

# ---------------- rational certificate: 2 real rational preimages ----------------
print("\n=== rational non-injectivity certificate ===")
# w=1, w=2, C=1:
Phi1 = Phi.subs(w, 1); Phi2 = Phi.subs(w, 2)
Bv = (Phi1 - Phi2)/(1 - 2)
Av = Bv*1 - Phi1
print(f"Phi(1) = {Phi1}, Phi(2) = {Phi2}  ->  target (A,B,C) = ({Av}, {Bv}, 1)")
def reconstruct(wstar, Bv, Av, Cv=sp.Rational(1)):
    gv = Bv*Cv - p.subs(w, wstar)          # gamma = BC - p(w)
    if gv == 0:
        return None
    xv = Cv/gv
    uv = wstar/gv
    yv = sp.simplify((uv - 1)/xv)
    tv = sp.simplify(gv - 1 - a_coef*(uv - 1))
    zv = sp.simplify(tv/xv**2)
    return (sp.simplify(xv), sp.nsimplify(yv), sp.nsimplify(zv))
pts = []
for wv in (sp.Rational(1), sp.Rational(2)):
    pt = reconstruct(wv, Bv, Av)
    im = [sp.simplify(f.subs({x: pt[0], y: pt[1], z: pt[2]})) for f in (f1, f2, f3)]
    print(f"  w* = {wv}: point {pt}  ->  {im}")
    pts.append(pt)
print("  distinct real points:", pts[0] != pts[1])
# remaining roots of the fiber quintic at this target
h = Phi - Bv*w + Av
rem = sp.factor(h/((w-1)*(w-2)))
print("  fiber quintic = (w-1)(w-2) *", rem)
rr = sp.nroots(sp.Poly(sp.expand(rem), w))
print("  remaining roots:", [complex(r2) for r2 in rr])

# ---------------- fold wall certificate: target (-1,-1,1) ----------------
print("\n=== fold certificate: (A,B,C) = (-1,-1,1), double root w=1 ===")
Bv2, Av2 = sp.Rational(-1), sp.Rational(-1)
hf = Phi - Bv2*w + Av2
print("  h_f(w) =", sp.factor(hf))
gam_at_double = Bv2*1 - p.subs(w, 1)
print("  gamma(w=1) = BC - p(1) =", gam_at_double, "  (double root escapes exactly)")
remf = sp.factor(hf/(w-1)**2)
rrf = sp.nroots(sp.Poly(sp.expand(remf), w))
print("  simple roots:", [complex(r2) for r2 in rrf])
for r2 in rrf:
    gv = Bv2 - p.subs(w, r2)
    xv = 1/gv; uv = r2/gv; yv = (uv-1)/xv; tv = gv - 1 - a_coef*(uv-1); zv = tv/xv**2
    s1, s2, s3 = [complex(f.subs({x: xv, y: yv, z: zv})) for f in (f1, f2, f3)]
    print(f"   w*={complex(r2):.4f}: preimage -> ({s1:.4f},{s2:.4f},{s3:.4f})")
