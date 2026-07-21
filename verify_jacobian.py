"""
Verification of the claimed counterexample to the Jacobian Conjecture
posted July 20, 2026 (L. Alpoge, crediting Claude Fable 5).

Claim: F: C^3 -> C^3 with
  f1 = (1+xy)^3 z + y^2 (1+xy) (4+3xy)
  f2 = y + 3x(1+xy)^2 z + 3 x y^2 (4+3xy)
  f3 = 2x - 3x^2 y - x^3 z
has (i) det JF = -2 (a nonzero constant) and (ii) is not injective,
with (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2) all mapping to (-1/4,0,0).
"""

import sympy as sp

x, y, z = sp.symbols("x y z")
s = 1 + x * y

f1 = s**3 * z + y**2 * s * (4 + 3 * x * y)
f2 = y + 3 * x * s**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
f3 = 2 * x - 3 * x**2 * y - x**3 * z

F = [f1, f2, f3]

# (i) Jacobian determinant
J = sp.Matrix(F).jacobian([x, y, z])
det = sp.expand(J.det())
print("det J(F) =", sp.factor(det))

# (ii) check the three alleged colliding points
pts = [
    (0, 0, sp.Rational(-1, 4)),
    (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
]
for p in pts:
    image = [sp.simplify(f.subs({x: p[0], y: p[1], z: p[2]})) for f in F]
    print("F%s =" % (p,), image)

# sanity: verify the map is genuinely polynomial and degrees
print("degrees:", [sp.Poly(f, x, y, z).total_degree() for f in F])
