"""
Stage 5c (fast): the Nagata-conjugated map G = sigma o F4 o sigma^{-1}, verified
POINT-WISE with exact rational arithmetic -- no giant symbolic expansion.

sigma (Nagata): w' = x z + y^2 ; (x - 2w'y - w'^2 z, y + w' z, z). det J sigma = 1,
sigma o sigma^{-1} = id (checked in stage 5).
Reports: degree bounds, det J G at 3 random rational points (exact), and the
non-injectivity transport: sigma(P1), sigma(P2) both map to sigma(target).
"""
import sympy as sp
from sympy import symbols, diff, integrate, expand, cancel, Matrix, Rational as R
import itertools

x, y, z = symbols("x y z")

def build(d):
    w = symbols("w")
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
    return f1, f2, f3

f1, f2, f3 = build(4)
F4 = [f1, f2, f3]
JF4 = [[sp.diff(f, v) for v in (x, y, z)] for f in F4]

w_p = x*z + y**2
sinv = [x + 2*w_p*y - w_p**2*z, y - w_p*z, z]
Jsinv = [[sp.diff(f, v) for v in (x, y, z)] for f in sinv]
w_s = x*z + y**2
sig = [x - 2*w_s*y - w_s**2*z, y + w_s*z, z]
Jsig = [[sp.diff(f, v) for v in (x, y, z)] for f in sig]

# --- degree bounds (deg sigma = 5): ---
print("degree of F4 components:", [sp.Poly(f, x, y, z).total_degree() for f in F4])
print("degree of sinv components:", [sp.Poly(f, x, y, z).total_degree() for f in sinv])
# F4 o sinv: deg <= 5*17 = 85, 80, 20 ; w' under G-composition: deg <= max(85+20, 2*80) = 160
# G1 <= max(85, 160+80=240, 2*160+20=340) = 340 ; G2 <= 160+20 = 180 ; G3 = 20
print("degree bounds for G = sigma o F4 o sinv: G1 <= 340, G2 <= 180, G3 = 20")

# --- det J G (pt) = det J sigma(F4(sinv pt)) * det J F4(sinv pt) * det J sinv(pt)
def ev(expr, pt):
    return sp.Rational(expr.subs({x: pt[0], y: pt[1], z: pt[2]}))

def detJ(J, pt):
    M = [[ev(J[i][j], pt) for j in range(3)] for i in range(3)]
    return sp.Matrix(M).det()

def apply_map(comp_fns, pt):
    return tuple(sp.simplify(ev(f, pt)) for f in comp_fns)

for pt in [(sp.Rational(1,3), sp.Rational(-2,5), sp.Rational(4,7)),
           (sp.Rational(2), sp.Rational(1,2), sp.Rational(-3)),
           (sp.Rational(-1,4), sp.Rational(3,2), sp.Rational(1,6))]:
    mid = apply_map(sinv, pt)
    top = apply_map(F4, mid)
    def subdict(pt):
        return {x: sp.Rational(pt[0]), y: sp.Rational(pt[1]), z: sp.Rational(pt[2])}
    Ms = sp.Matrix([[sp.Rational(Jsig[i][j].subs(subdict(top))) for j in range(3)] for i in range(3)])
    Mf = sp.Matrix([[sp.Rational(JF4[i][j].subs(subdict(mid))) for j in range(3)] for i in range(3)])
    Mi = sp.Matrix([[sp.Rational(Jsinv[i][j].subs(subdict(pt))) for j in range(3)] for i in range(3)])
    det_prod = Ms.det() * Mf.det() * Mi.det()
    print(f"  pt={pt}: detJG = detJsigma(F4s)*detJF4(sinv)*detJsinv = "
          f"{Ms.det()} * {Mf.det()} * {Mi.det()} = {det_prod}")

# --- non-injectivity transport ---
P1 = (R(-5,26), R(31,5), R(-714922,3375))
P2 = (R(5,46), R(-36,5), R(226228,375))
target = (-sp.Rational(31,5),)*2 + (sp.Rational(1),)
def subdict(pt):
    return {x: sp.Rational(pt[0]), y: sp.Rational(pt[1]), z: sp.Rational(pt[2])}
def apply_exact(fns, pt):
    return tuple(sp.Rational(f.subs(subdict(pt))) for f in fns)
sP1 = apply_exact(sig, P1); sP2 = apply_exact(sig, P2)
GsP1 = apply_exact(sig, apply_exact(F4, P1))
GsP2 = apply_exact(sig, apply_exact(F4, P2))
print("  sigma(F4(P1)) =", GsP1)
print("  sigma(F4(P2)) =", GsP2, " equal:", GsP1 == GsP2, " distinct inputs:", sP1 != sP2)
print("  preimages of sigma(target) under G = sig o F4 o sinv should be sig(P_i):")
for Pi in (P1, P2):
    qi = apply_exact(sig, Pi)                      # = sigma(P_i)
    Gi = apply_exact(sig, apply_exact(F4, apply_exact(sinv, qi)))
    print("   G(sigma(P_i)) =", Gi, " == sigma(target):",
          Gi == apply_exact(sig, target))
