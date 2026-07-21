"""
Deep-dive verification of the Alpoge/Fable Jacobian counterexample.
Part 1: structure, hidden symmetries, the escape-to-infinity curve,
and the 'fault divisor'.

Posted map F: C^3 -> C^3:
  f1 = (1+xy)^3 z + y^2 (1+xy) (4+3xy)
  f2 = y + 3x(1+xy)^2 z + 3 x y^2 (4+3xy)
  f3 = 2x - 3x^2 y - x^3 z
"""
import sympy as sp

x, y, z, t = sp.symbols("x y z t")

# ---- A factored normal form ---------------------------------------------
# The posted polynomials factor beautifully:
#   with w = 1+xy, A = w^2 z + y^2 (4+3xy), B = 3y + xz:
#   f1 = w*A, f2 = y + 3x*A, f3 = 2x - x^2 B
w = 1 + x * y
A = w**2 * z + y**2 * (4 + 3 * x * y)
B = 3 * y + x * z
f1 = sp.expand(w * A)
f2 = sp.expand(y + 3 * x * A)
f3 = sp.expand(2 * x - x**2 * B)

posted1 = sp.expand(w**3 * z + y**2 * w * (4 + 3 * x * y))
posted2 = sp.expand(y + 3 * x * w**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
posted3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
print("factored form matches posting:",
      sp.expand(f1 - posted1) == 0 and sp.expand(f2 - posted2) == 0
      and sp.expand(f3 - posted3) == 0)

# ---- The Jacobian test ----------------------------------------------------
F = sp.Matrix([f1, f2, f3])
det = sp.factor(F.jacobian([x, y, z]).det())
print("det J(F) =", det)

# ---- Hidden odd symmetry: F(-x,-y,z) = (f1, -f2, -f3) ----------------------
flip = {x: -x, y: -y}
sym = [sp.expand(f1.subs(flip) - f1) == 0,
       sp.expand(f2.subs(flip) + f2) == 0,
       sp.expand(f3.subs(flip) + f3) == 0]
print("odd symmetry  F(-x,-y,z) = (f1,-f2,-f3):", all(sym))

# ---- The escape curve: properness fails ------------------------------------
esc = {x: t, y: -1 / t, z: 5 / t**2}
print("F(t, -1/t, 5/t^2) =", [sp.simplify(f.subs(esc)) for f in (f1, f2, f3)],
      "   <-- point escapes to infinity, image tends to (0,0,0)")

# ---- The fault divisor xy = -1 ---------------------------------------------
on_div = {y: -1 / x}
print("F(x, -1/x, z)  =", [sp.simplify(f.subs(on_div)) for f in (f1, f2, f3)],
      "   <-- divisor {xy=-1} maps bijectively onto {0} x C* x C")

# ---- Fiber over the origin: exactly one point (symbolic argument) ----------
# f1 = w*A, f2 = y + 3x*A, so f1 = f2 = 0 gives two cases:
#   case w = 0:  then A = y^2 (since 4+3xy = 1) and f2 = y + 3x y^2 = -2y != 0. Contradiction.
Aw0 = sp.expand(A.subs(y, -1 / x))
f2w0 = sp.expand((y + 3 * x * A).subs({y: -1 / x}))
print("on w=0:  A =", Aw0, ", f2 =", f2w0, " (nonzero for x in C*, so no zeros on the divisor)")
#   case A = 0:  then f2 = y = 0, hence w = 1, A = z, so z = 0, and f3 = 2x = 0.
Aw_y0 = sp.expand(A.subs(y, 0))
print("on A=0, y=0: w =", sp.expand(w.subs(y, 0)), ", A =", Aw_y0,
      ", f3 =", sp.expand(f3.subs({y: 0, z: 0})))

# ---- The three known collision points ---------------------------------------
p1 = (0, 0, sp.Rational(-1, 4))
p2 = (1, sp.Rational(-3, 2), sp.Rational(13, 2))
p3 = (-1, sp.Rational(3, 2), sp.Rational(13, 2))
for p in (p1, p2, p3):
    print(f"F{p} =", [sp.simplify(f.subs({x: p[0], y: p[1], z: p[2]})) for f in (f1, f2, f3)])
print("note: p2 and p3 are related by the odd symmetry (flip x,y), p1 is a fixed point of it")
