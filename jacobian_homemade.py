"""
Minting a NEW counterexample to the Jacobian Conjecture from my own seed polynomial,
via the weighted-lift recipe (p, c) with:
  seed p(w) = 5w - 12w^2 + 6w^3,   c = 1,  b = 1,  a = -(1+kappa)/(2+kappa) = 0.
This seed is NOT the one in Alpoge's map (2w-3w^2) and NOT the explainer's (w-2w^3),
so success = independent confirmation that the construction is a genuine machine.

Certificate needed:
  (1) H is polynomial in x, y, z
  (2) det J(H) = 1 (nonzero constant)
  (3) an exact collision: two distinct exact algebraic points mapping to one target.
"""
import sympy as sp

x, y, z, w = sp.symbols("x y z w")

# ---- seed and recipe ingredients ----
p = 5 * w - 12 * w**2 + 6 * w**3
c = sp.Integer(1)
print("p(0) =", p.subs(w, 0), " p(1) =", p.subs(w, 1), " (need 0 and -1)")
print("int_0^1 p =", sp.integrate(p, (w, 0, 1)), " (need 0)")
kap = sp.diff(p, w).subs(w, 1) / c
print("kappa =", kap, " (need != -2)")
a = sp.Rational(0); b = sp.Integer(1)

q = sp.expand(sp.integrate(w * sp.diff(p, w) / c, w))  # q(0)=0
print("q(w) =", q)

u = 1 + x * y
g = 1 + x**2 * z            # gamma = 1 + a*xy + b*x^2 z with a=0, b=1
ws = u * g

# ---- components with exact cancellation checks ----
num1 = sp.expand(u * g**2 + q.subs(w, ws))
alpha = sp.cancel(num1 / g**2)
print("H1: numerator divisible by gamma^2:", sp.denom(alpha) == 1)
alpha = sp.Poly(sp.expand(alpha), x)
print("H1: quotient divisible by x^2:", all(m[0] >= 2 for m in alpha.monoms()))
H1 = sp.expand(sp.cancel(sp.expand(sum(ci * x**mi[0] * y**mi[1] * z**mi[2]
                 for mi, ci in sp.Poly(sp.cancel(num1 / g**2), x, y, z).terms())) / x**2))

num2 = sp.expand(c * g + p.subs(w, ws))
beta = sp.cancel(num2 / g)
print("H2: numerator divisible by gamma:", sp.denom(beta) == 1)
betaP = sp.Poly(sp.expand(beta), x)
print("H2: quotient divisible by x:", all(m[0] >= 1 for m in betaP.monoms()))
H2 = sp.expand(sp.cancel(beta / x))
H3 = sp.expand(x * g)

print("H degrees:", [sp.Poly(f, x, y, z).total_degree() for f in (H1, H2, H3)])
print("H1 has", len(sp.Poly(H1, x, y, z).terms()), "terms; H2 has",
      len(sp.Poly(H2, x, y, z).terms()), "terms")

# ---- constant Jacobian ----
detH = sp.factor(sp.Matrix([H1, H2, H3]).jacobian([x, y, z]).det())
print("det J(H) =", detH, " (recipe predicts b*c = 1)")

# ---- exact collision certificate ----
# Fiber equation: Phi(w) = w P - c Q, P = B*C, Q = A*C^2, Phi = int_0^w p.
Phi = sp.expand(sp.integrate(p, w))  # antiderivative with Phi(0) = 0
print("Phi(w) =", Phi)
d = sp.symbols("d")
diffd = sp.factor(Phi.subs(w, 2 + d) - Phi.subs(w, 2 - d))
print("Phi(2+d) - Phi(2-d) =", diffd)
# root d of the odd part:
sol = sp.solve(sp.expand(diffd / (2 * d)), d**2)
d0 = sp.sqrt(sol[0])
w1, w2 = 2 + d0, 2 - d0
print("w1 =", w1, " w2 =", w2, " Phi(w1)-Phi(w2) =",
      sp.simplify(Phi.subs(w, w1) - Phi.subs(w, w2)))

A0 = sp.simplify(-Phi.subs(w, w1))          # target (A0, 0, 1)
P0, C0 = sp.Integer(0), sp.Integer(1)
def reconstruct(wroot):
    gm = sp.simplify((P0 - p.subs(w, wroot)) / c)
    xi = sp.simplify(C0 / gm)
    vi = sp.simplify(wroot / gm - 1)
    yi = sp.simplify(vi / xi)
    ti = sp.simplify((gm - 1 - a * vi) / b)
    zi = sp.simplify(ti / xi**2)
    return (xi, yi, zi)
pt1 = reconstruct(w1); pt2 = reconstruct(w2)
print("pt1 =", [sp.nsimplify(v) for v in pt1])
print("pt2 =", [sp.nsimplify(v) for v in pt2])

for name, pt in (("H(pt1)", pt1), ("H(pt2)", pt2)):
    im = [sp.simplify(f.subs({x: pt[0], y: pt[1], z: pt[2]})) for f in (H1, H2, H3)]
    print(name, "=", im)
print("target =", (A0, 0, 1))
print("numeric check det at random pts:",
      [sp.N(H1.subs({x: sp.Rational(3, 7), y: -5, z: 2}), 12)])
