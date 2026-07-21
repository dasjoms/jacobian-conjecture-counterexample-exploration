"""
Note 10, stage F: first-ever determinant check for F8 (nobody computed det JF8 before -
d=6,7 chambers relied on the recipe). Pointwise-EXACT: build the recipe's rational
expression tree (no expansion), differentiate the tree, evaluate at rational points.
Recipe predicts det JF8 = b*c = 1 identically.
"""
import sympy as sp, random
x, y, z, w = sp.symbols("x y z w")
d = 8
p = 2*w - 3*w**2 + w*(1-w)*(w**(d-2) - sp.Rational(6, d*(d+1)))
p = sp.expand(p)
c = sp.Integer(1); b = sp.Integer(1)
q = sp.expand(sp.integrate(w*sp.diff(p, w)/c, w))
kap = sp.Rational(sp.diff(p, w).subs(w, 1))/c
a = sp.Rational(-(1+kap)/(2+kap))
u = 1 + x*y
g = 1 + a*x*y + b*x**2*z
ws = u*g
alpha = u + q.subs(w, ws)/g**2
beta = c + p.subs(w, ws)/g
f1 = alpha/x**2; f2 = beta/x; f3 = x*g
print("recipe built: kappa =", kap, " a =", a, flush=True)
J = sp.Matrix([f1, f2, f3]).jacobian([x, y, z])   # tree differentiation (fast)
print("jacobian tree done", flush=True)
random.seed(8)
pts = []
for _ in range(5):
    while True:
        pt = {x: sp.Rational(random.randint(-4,4), random.randint(1,5)),
              y: sp.Rational(random.randint(-4,4), random.randint(1,5)),
              z: sp.Rational(random.randint(-4,4), random.randint(1,5))}
        gv = sp.simplify(g.subs(pt))
        if pt[x] != 0 and gv != 0: break
    pts.append(pt)
ok = 0
for i, pt in enumerate(pts):
    vals = [sp.nsimplify(entry.subs(pt), rational=False) for entry in J]  # exact rationals
    det = sp.Matrix(3, 3, lambda ii, jj: vals[3*ii+jj]).det()
    print(f"  point {i}: det = {sp.factor(det)}", flush=True)
    if sp.simplify(det - 1) == 0: ok += 1
print(f"det == 1 at {ok}/5 rational points (exact arithmetic)")
