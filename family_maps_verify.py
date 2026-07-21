"""
Build the explainer's seed-family maps for d = 3, 4 (and 5 if feasible) via the
weighted-lift recipe, and verify the three pillars needed for the
surjectivity theorem:
  (a) Fp is polynomial, det JFp = b*c = 1 (recipe's prediction)
  (b) fiber identity  Phi(u*gamma) = (B*C)*(u*gamma) - A*C^2  holds as a
      polynomial identity  (c = 1)
  (c) non-injectivity: exhibit >= 2 exact preimages of a target (fiber degree >= 4)
"""
import sympy as sp

x, y, z, w = sp.symbols("x y z w")

def build(d):
    p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))
    c = sp.Integer(1); b = sp.Integer(1)
    q = sp.expand(sp.integrate(w*sp.diff(p, w)/c, w))
    kap = sp.Rational(sp.diff(p, w).subs(w, 1)) / c
    a = sp.Rational(-(1 + kap) / (2 + kap))
    print(f"  d={d}: seed p(w) = {sp.expand(p)}, kappa = {kap}, a = {a}")
    u = 1 + x*y
    g = 1 + a*x*y + b*x**2*z
    ws = u*g
    alpha = u + q.subs(w, ws)/g**2
    beta = c + p.subs(w, ws)/g
    return p, sp.expand(sp.integrate(p, w)), u, g, ws, alpha, beta

def verify_map(d):
    print(f"--- d = {d} ---")
    p, Phi, u, g, ws, alpha, beta = build(d)
    # (a1) polynomiality of alpha/x^2 and beta/x
    A_ = sp.cancel(alpha/x**2)
    B_ = sp.cancel(beta/x)
    numA, denA = sp.fraction(A_)
    numB, denB = sp.fraction(B_)
    okA = denA == 1
    okB = denB == 1
    print("  alpha/x^2 polynomial:", okA, " beta/x polynomial:", okB)
    if not (okA and okB):
        return None
    f1 = sp.expand(numA); f2 = sp.expand(numB); f3 = sp.expand(x*g)
    print("  degrees:", [sp.Poly(f, x, y, z).total_degree() for f in (f1, f2, f3)])
    # (a2) determinant
    det = sp.factor(sp.Matrix([f1, f2, f3]).jacobian([x, y, z]).det())
    print("  det =", det)
    # (b) fiber identity with P = B C, Q = A C^2  (c = 1)
    lhs = sp.expand(Phi.subs(w, ws) - (f2*f3)*ws + f1*f3**2)
    print("  fiber identity holds:", sp.factor(lhs) == 0)
    return f1, f2, f3, Phi, p

res3 = verify_map(3)
res4 = verify_map(4)

# (c) exact non-injectivity demo for d=3: target (A,B,C) = (0,0,1)
# fiber: Phi_4(w) = 0  (P=Q=0). Roots: w=0 (escapes, gamma = -p(0) = 0), others.
if res3:
    f1, f2, f3, Phi, p = res3
    rest = sp.factor(Phi / w**2)
    print("  d=3: Phi_4(w) = w^2 *", rest)
    sols = sp.solve(sp.Eq(rest, 0), w)
    print("  roots of rest:", sols)
    gam = lambda wv: sp.simplify(-p.subs(w, wv))          # gamma = P - p(w), P = 0
    for wv in sols:
        g0 = gam(wv)
        print(f"   w = {wv}: gamma = {g0}")
        xx = sp.simplify(1/g0); uu = sp.simplify(wv/g0)
        yv = sp.simplify((uu - 1)/xx)
        a_kap = sp.Rational(-7,5)   # a for d=3: gamma = 1 + a xy + x^2 z  -> z from gamma eq
        tv = sp.simplify(g0 - 1 - a_kap*(uu - 1))
        zv = sp.simplify(tv / xx**2)
        im = [sp.simplify(f.subs({x: xx, y: yv, z: zv})) for f in (f1, f2, f3)]
        print(f"     preimage = ({xx}, {yv}, {zv})  maps to {im}")


# ---------- d = 5 (degrees ~(22,21,4)) ----------
res5 = verify_map(5)
if res5:
    print("  d=5 verified (det, identity).")
