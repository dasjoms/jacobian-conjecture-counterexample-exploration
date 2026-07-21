"""
Note 6, stage A: the wall discriminant for F4's fiber family, the F sanity check,
and the special strata (cusps from p' roots; bitangent nodes).

Fiber family: h(w; s, r) = Phi(w) - s w + r,  (s = BC, r = AC^2).
Wall:  Disc(s, r) = resultant(h, h', w) = 0.
F family (sanity): h_F = -w^3 + w^2 - P w + Q ; expect 4P^3 - P^2 - 18 P Q + 27 Q^2 + 4Q
(note-2 wall), possibly up to a nonzero scalar.
"""
import sympy as sp
from sympy import symbols, diff, integrate, expand, factor, resultant, Rational as R

w, s, r, P_, Q_ = symbols("w s r P Q")

p4 = -w**4 + w**3 - R(27,10)*w**2 + R(17,10)*w
Phi4 = expand(integrate(p4, w))

h4 = Phi4 - s*w + r
print("computing Disc_4(s,r) = resultant(h, h', w) ...", flush=True)
D4 = sp.expand(resultant(h4, diff(h4, w), w))
D4 = sp.factor(D4 / sp.factor(D4).args[0]) if False else D4
print("degree in (s,r):", sp.Poly(D4, s, r).total_degree(), " terms:", len(sp.Poly(D4, s, r).terms()))
print("content:", sp.factor(sp.gcd_list([c for c in sp.Poly(D4, s, r).coeffs()])))
D4x = sp.expand(D4 / sp.gcd_list(list(sp.Poly(D4, s, r).coeffs())))
print("Disc_4(s,r) =", D4x)
print("Disc_4(0,0) =", D4x.subs({s: 0, r: 0}), "  -> C=0 plane lies in wall, as predicted")
print("irreducible factorization:", sp.factor(D4x) == D4x)

# --- sanity: F's cubic family ---
hF = -w**3 + w**2 - P_*w + Q_
DF = sp.expand(resultant(hF, diff(hF, w), w))
print("\nF-wall Disc(P,Q) =", sp.factor(DF))
target = 4*P_**3 - P_**2 - 18*P_*Q_ + 27*Q_**2 + 4*Q_
print("matches note-2 R(P,Q) up to scalar:", sp.factor(DF) , "/", target, "ratio:",
      sp.simplify(DF/target))

# --- cusp strata: roots of p4' (tangent lines of order >= 3) ---
print("\n--- (3)-cusp points of the wall ---")
roots_pp = sp.nroots(diff(p4, w), n=40)
for wv in roots_pp:
    sv = p4.subs(w, wv); rv = wv*sv - Phi4.subs(w, wv)
    tag = "REAL" if abs(sp.im(wv)) < 1e-20 else "complex"
    print(f"  w0={complex(wv):.8f} -> (s,r)=({complex(sv):.8f},{complex(rv):.8f})  [{tag}]")

# --- bitangent nodes: solve p(w1)=p(w2), Phi-secant slope = p(w1), off diagonal ---
print("\n--- (2,2) bitangent nodes ---")
w1, w2 = symbols("w1 w2")
eq1 = sp.cancel((p4.subs(w, w2) - p4.subs(w, w1))/(w2 - w1))
eq2 = sp.cancel((Phi4.subs(w, w2) - Phi4.subs(w, w1))/(w2 - w1) - p4.subs(w, w1))
eq1 = sp.expand(eq1); eq2 = sp.expand(eq2)
print("  eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree())
print("  computing lex GB over Q[w1,w2] ...", flush=True)
gb = sp.groebner([eq1, eq2], w1, w2, order="lex")
print("  GB basis:", [sp.Poly(g, w1, w2).total_degree() for g in gb.polys], "polynomials")
sols = []
if gb.polys:
    univ = gb.polys[0].expr
    for wv in sp.nroots(sp.Poly(univ, w1), n=30, maxsteps=300):
        for g in gb.polys[1:]:
            pass
    print("  univariate eliminant degree:", sp.Poly(univ, w1).degree())
    # solve triangular via sympy.solve is awkward; use nsolve from backsolving
print("done stage A part 1")
