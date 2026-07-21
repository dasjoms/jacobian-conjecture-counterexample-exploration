"""
Stage 5: loose ends.
(a) third REAL preimage at (-31/5,-31/5,1) (w* = the real cubic root)
(b) whisker check: real targets near the missing curve M_R are attained
(c) Nagata conjugation of F4: automorphism sanity, degree/sparsity table,
    det = 1 via exact rational point evaluation; non-injectivity transport
"""
import sympy as sp
from sympy import symbols, diff, integrate, expand, cancel, Matrix, Rational as R, nsimplify
import math

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
    return sp.expand(p), sp.expand(sp.integrate(p, w)), a, f1, f2, f3

p, Phi, a_coef, f1, f2, f3 = build(4)

def reconstruct(wstar, Bv, Av, Cv=R(1)):
    gv = Bv*Cv - p.subs(w, wstar)
    if gv == 0: return None
    xv = Cv/gv; uv = wstar/gv
    yv = (uv - 1)/xv
    zv = (gv - 1 - a_coef*(uv - 1))/xv**2
    return (xv, yv, zv)

print("=== (a) third real preimage at (-31/5,-31/5,1) ===")
Bv, Av = R(-31,5), R(-31,5)
h = Phi - Bv*w + Av
rem = sp.Poly(sp.quo(h, (w-1)*(w-2), w), w)
rts = sp.nroots(rem, n=40)
w3 = [r for r in rts if abs(sp.im(r)) < 1e-25][0]
pt = reconstruct(sp.nsimplify(w3), Bv, Av)
im = [complex(f.subs({x: pt[0], y: pt[1], z: pt[2]})) for f in (f1, f2, f3)]
print(f"  w* = {complex(w3)} (real)  gamma = {Bv - p.subs(w, w3)}")
print(f"  preimage ~ ({[complex(v) for v in pt]})  -> {im}")

print("\n=== (b) whisker: targets near M_R are attained over R ===")
crit = [r for r in sp.nroots(diff(p, w), n=40) if abs(sp.im(r)) < 1e-25][0]
s0 = p.subs(w, crit); r0 = crit*s0 - Phi.subs(w, crit)
print(f"  cusp target M_R point: A={float(r0):.8f} B={float(s0):.8f} C=1  (MISSED over R)")
for eps, dA, dB in [(1e-3,1,0), (1e-3,-1,0), (1e-3,0,1), (1e-3,0,-1), (1e-2,1,1)]:
    Av2 = r0 + dA*eps; Bv2 = s0 + dB*eps
    hh = phi_h = sp.expand(Phi - Bv2*w + Av2)
    roots = sp.nroots(sp.Poly(hh.evalf(40), w), n=40, maxsteps=400)
    real_simple = []
    for rr in roots:
        if abs(sp.im(rr)) < 1e-12:
            # simplicity: |h'(rr)|
            dv = abs(diff(p, w).subs(w, rr) - Bv2)
            real_simple.append(float(dv) > 1e-8)
    attained = any(real_simple)
    print(f"   dA={dA*eps:+.3f} dB={dB*eps:+.3f}: real roots={sum(1 for rr in roots if abs(sp.im(rr))<1e-12)},"
          f" real&simple={sum(real_simple)}  -> over-R attained: {attained}")

print("\n=== (c) Nagata orbit ===")
X, Y, Z = symbols("X Y Z")
wxyz = X*Z + Y**2
sx = X - 2*wxyz*Y - wxyz**2*Z
sy = Y + wxyz*Z
sz = Z
wxz2 = x*z + y**2
sinv = [x + 2*(x*z+y**2)*y - (x*z+y**2)**2*z, y - (x*z+y**2)*z, z]
sig = [sx, sy, sz]
# check sigma o sinv = id
comp = [sp.expand(si.subs({X: sinv[0], Y: sinv[1], Z: sinv[2]})) for si in sig]
print("  sigma o sinv =", comp, " (identity: ", comp == [x, y, z], ")")
Jsig = Matrix([sx, sy, sz]).jacobian([X, Y, Z]).det()
print("  det J sigma =", sp.factor(Jsig))

print("  composing G = sigma o F4 o sinv ...", flush=True)
sub = {x: sinv[0], y: sinv[1], z: sinv[2]}
F4s = [sp.expand(f.subs(sub)) for f in (f1, f2, f3)]
print("   F4 o sinv degrees:", [sp.Poly(f, x, y, z).total_degree() for f in F4s],
      "terms:", [len(sp.Poly(f, x, y, z).terms()) for f in F4s])
G = [sp.expand(sig[i].subs({X: F4s[0], Y: F4s[1], Z: F4s[2]})) for i in range(3)]
print("   G degrees:", [sp.Poly(f, x, y, z).total_degree() for f in G],
      "terms:", [len(sp.Poly(f, x, y, z).terms()) for f in G])
