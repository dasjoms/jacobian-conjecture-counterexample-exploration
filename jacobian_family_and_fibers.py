"""
Independent verification of the CLAIMED second counterexample (degree-4 fiber map G),
numeric fiber census for both maps, and mod-p bijection tests.

G (from jacobianfun.org, unverified until checked here):
  u = 1 + 3xy, gamma = 1 - 4xy - x^2 z
  G = ( (2u + u^2 - 3 u^4 gamma^2)/x^2,  (1 + u - 2 u^3 gamma^2)/x,  x gamma )
Claims: (a) G is polynomial, (b) det JG = -6, (c) G(1,0,0) = G(-1,0,2) = (0,0,1).
"""
import sympy as sp
import numpy as np

x, y, z = sp.symbols("x y z")

# ================= verify G =================
u = 1 + 3 * x * y
g = 1 - 4 * x * y - x**2 * z
N1 = sp.expand(2 * u + u**2 - 3 * u**4 * g**2)
N2 = sp.expand(1 + u - 2 * u**3 * g**2)

# (a) polynomiality: no x^0 or x^1 terms in N1, no x^0 terms in N2
ok1 = all(m[0] >= 2 for m in sp.Poly(N1, x).monoms())
ok2 = all(m[0] >= 1 for m in sp.Poly(N2, x).monoms())
print("G1 numerator divisible by x^2:", ok1, " G2 numerator divisible by x:", ok2)
g1 = sp.expand(sp.cancel(N1 / x**2))
g2 = sp.expand(sp.cancel(N2 / x))
g3 = sp.expand(x * g)
print("G component degrees:", [sp.Poly(f, x, y, z).total_degree() for f in (g1, g2, g3)])

# (b) constant Jacobian
detG = sp.factor(sp.Matrix([g1, g2, g3]).jacobian([x, y, z]).det())
print("det J(G) =", detG)

# (c) collision certificate
for pt in [(1, 0, 0), (-1, 0, 2)]:
    print(f"G{pt} =", [sp.Integer(f.subs({x: pt[0], y: pt[1], z: pt[2]})) for f in (g1, g2, g3)])

# ============= exact 3-point fiber of the original map at (-8/27, 0, 1) ======
# per explainer: (-3/2, 4/3, 104/27) and a conjugate radical pair
w0 = 1 + x * y
A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
f1 = sp.expand(w0 * A0); f2 = sp.expand(y + 3 * x * A0); f3 = sp.expand(2 * x - x**2 * (3 * y + x * z))
r = sp.sqrt(3)
pts = [(-sp.Rational(3, 2), sp.Rational(4, 3), sp.Rational(104, 27)),
       (sp.Rational(3, 4), 2 * (-1 + r) / 3, 8 * (13 - 9 * r) / 27),
       (sp.Rational(3, 4), -2 * (1 + r) / 3, 8 * (13 + 9 * r) / 27)]
for pt in pts:
    im = [sp.simplify(f.subs({x: pt[0], y: pt[1], z: pt[2]})) for f in (f1, f2, f3)]
    print(f"F{pt} = {im}  (all should equal (-8/27, 0, 1))")
