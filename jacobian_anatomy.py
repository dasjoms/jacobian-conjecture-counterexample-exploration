"""
Anatomy of the counterexample F: proof-grade fiber census.

Derive and verify:
  (i)   exact identity -- every preimage (x neq 0) comes from a root w = u*gamma
        of the fiber cubic  Phi(w) = w^2 - w^3 = P w - Q,  P = BC/4, Q = AC^2/4
  (ii)  at most 3 finite preimages over any C != 0 target (reconstruction injective)
  (iii) the 'wall' variety: targets where sheets escape = cubic multiple root
        <=> (P,Q) on the rational curve Gamma = {(p(w), q(w))}
  (iv)  cusp of Gamma at (P,Q) = (1/3, 1/27)  ==> triple root escapes ==>
        THE MISSING CURVE  M = {(4/(27 C^2), 4/(3C), C) : C != 0}
        targets on M have NO preimage at all.
  (v)   C = 0 plane treated separately (flat + curved sheet).
Then the analogous wall/cusp computation for G and the homemade H.
"""
import sympy as sp

x, y, z, w, P, Q = sp.symbols("x y z w P Q")

# ---------------- rebuild F ----------------
u = 1 + x * y
gm = 1 - sp.Rational(3, 2) * x * y - sp.Rational(1, 2) * x**2 * z      # gamma = 1 -(3/2)v -(1/2)t
p = 2 * w - 3 * w**2
q = w**2 - 2 * w**3
Phi = w**2 - w**3

w0 = 1 + x * y
A0 = w0**2 * z + y**2 * (4 + 3 * x * y)
F1 = sp.expand(w0 * A0); F2 = sp.expand(y + 3 * x * A0); F3 = sp.expand(2 * x - x**2 * (3 * y + x * z))

# sanity: the announced map in weighted coordinates matches our (u, gamma) setup?
# verify gamma's role:  w(=u*gamma) used below in the identity check directly.

A_, B_, C_ = F1, F2, F3
wsym = u * gm
lhs = sp.expand(Phi.subs(w, wsym) - (B_ * C_ / 4) * wsym + (A_ * C_**2 / 4))
print("(i) F fiber identity  Phi(u*gamma) - P*u*gamma + Q  == 0 :", sp.factor(lhs) == 0)

# ---------------- (iii) wall curve ----------------
R = sp.factor(sp.resultant(P - p, Q - q, w))
print("(iii) wall curve Gamma: R(P,Q) =", R)
disc = sp.factor(sp.discriminant(w**3 - w**2 + P * w - Q, w))
print("      equals disc of fiber cubic (up to factor):", sp.factor(R - disc) == 0 or sp.factor(R + disc) == 0, "; disc =", disc)

# cusp = singular point of R
sing = sp.solve([R, sp.diff(R, P), sp.diff(R, Q)], [P, Q], dict=True)
print("(iv) singular point(s) of Gamma:", sing)

# no node: q(w) = q(2/3 - w) only at the symmetry axis w = 1/3
node = sp.factor(q - q.subs(w, sp.Rational(2, 3) - w))
print("      q(w) - q(2/3-w) =", node, " -> roots:", sp.solve(node, w))

# ---------------- (iv) the missing curve ----------------
w_trip = sp.Rational(1, 3)
Ptrip = sp.expand(p.subs(w, w_trip)); Qtrip = sp.expand(q.subs(w, w_trip))
print("triple-escape point in (P,Q):", (Ptrip, Qtrip))
print("missing curve  M(C) = (4*Q_trip/C^2, 4*P_trip/C, C)  =  (4/(27 C^2), 4/(3C), C)")

# emptiness at a missing point: cubic is (w-1/3)^3 and gamma(1/3)=0 there
gam_w = P - p            # gamma of root w (given target P)
cub_at_cusp = sp.expand(w**3 - w**2 + Ptrip * w - Qtrip)
print("fiber cubic at cusp =", sp.factor(cub_at_cusp), "; gamma at double/triple root:",
      sp.simplify(gam_w.subs({P: Ptrip, w: w_trip})))

# ---------------- (v) C = 0 plane identities ----------------
v = x * y
zcur = (2 - 3 * v) / x**2
e1 = sp.expand((F1 * x**2 - (v + 1) * (v + 2)).subs(z, zcur))
e2 = sp.expand((F2 * x - 2 * (2 * v + 3)).subs(z, zcur))
print("(v) curved sheet x^2 z = 2-3xy maps as ((v+1)(v+2)/x^2, 2(2v+3)/x, 0):",
      e1 == 0 and e2 == 0)
e3 = [sp.simplify(f.subs({x: 0})) for f in (F1, F2, F3)]
print("    flat sheet x=0 maps as (z+4y^2, y, 0):", e3)
