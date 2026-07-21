"""
Note 10, stage 1: (A) the d=3 chamber anchor (fiber 4) + (B) generic-seed harness.
p_d family re-derived; generic normalized deg-4 seed pipeline:
  seed p (deg 4, p(1) = -1, Phi(1) = 0) -> wall D(s,r) -> cusps/nodes/budget.
"""
import sympy as sp, json
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd
import numpy as np

w, s, r = symbols("w s r")

def wall_of(p):
    """exact primitive wall discriminant for a normalized seed p (int-scaled)."""
    Phi = expand(integrate(p, w))
    h = sp.expand(sp.together((Phi - s*w + r)*sp.factor(sp.denom(sp.together(Phi - s*w + r)))))
    # clear denominators robustly:
    coeffs = sp.Poly(expand(Phi - s*w + r), w).all_coeffs()
    L = 1
    for c in coeffs:
        L = sp.ilcm(L, sp.denom(c))
    h = expand(L*(Phi - s*w + r))
    D = resultant(h, diff(h, w), w)
    cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
    return sp.expand(D/cont), Phi

def cusps_of(p, Phi):
    out = []
    for wv in sp.nroots(diff(p, w), n=50):
        sv = p.subs(w, wv); rv = wv*sv - Phi.subs(w, wv)
        out.append((complex(sp.N(wv,50)), complex(sp.N(sv,50)), complex(sp.N(rv,50))))
    return out

def nodes_of(p, Phi, prec=60):
    """bitangent pairs via contact-image clustering (note-9 scheme, GB-free)."""
    w1, w2 = symbols("w1 w2")
    tt = symbols("tt")
    pa = p.subs(w, tt); Ph = Phi.subs(w, tt)
    # eliminate t from image equations s = p(t), r = t p(t) - Phi(t):
    # nodes = pairs with same image. Sample contacts numerically is not enough;
    # instead: resultant of (p(t1)-p(t2))/(t1-t2) & secant-eq, then cluster images.
    eq1 = sp.expand(sp.cancel((p.subs(w, w2) - p.subs(w, w1))/(w2 - w1)))
    eq2 = sp.expand(sp.cancel((Phi.subs(w, w2) - Phi.subs(w, w1))/(w2 - w1) - p.subs(w, w1)))
    gb = sp.groebner([eq1, eq2], w1, w2, order="lex")
    elim = None
    for g in gb.polys:
        if set(g.free_symbols) == {w2}:
            elim = g.expr
    fel = sp.factor(elim)
    sqpart = sp.sqf_part(sp.Poly(elim, w2).as_expr()) if False else None
    # use the SQUAREFREE part: kills the p'^2 diagonal content -> bitangent contacts
    cs, elimsf = sp.Poly(elim, w2).sqf_list()
    bit_elim = None
    for fac, mult in elimsf:
        if mult == 1 and sp.degree(fac) > 0:
            # expect the unique squarefree factor of degree (n-2)(n-3)
            if bit_elim is None or sp.degree(fac) > sp.degree(bit_elim):
                bit_elim = fac
    roots = sp.nroots(sp.Poly(bit_elim, w2), n=prec, maxsteps=3000)
    imgs = []
    for cv in roots:
        sv = complex(sp.N(p.subs(w, cv), prec)); rv = complex(sp.N(cv*p.subs(w, cv) - Phi.subs(w, cv), prec))
        imgs.append((complex(sp.N(cv, prec)), sv, rv))
    pairs, used, gaps = [], set(), []
    for i in range(len(imgs)):
        if i in used: continue
        best, bj = 1e30, -1
        for j in range(i+1, len(imgs)):
            if j in used: continue
            d = abs(imgs[i][1]-imgs[j][1]) + abs(imgs[i][2]-imgs[j][2])
            if d < best: best, bj = d, j
        used.add(i); used.add(bj)
        pairs.append((imgs[i], imgs[bj])); gaps.append(best)
    return pairs, gaps, elim

print("=== (A) d=3 chamber anchor (fiber 4) ===")
p3 = 2*w - 3*w**2 + w*(1-w)*(w - R(1,2))
p3 = sp.expand(p3)
print("p3 =", p3, " p3(1) =", p3.subs(w,1))
D3, Phi3 = wall_of(p3)
print("Phi3 =", Phi3, " Phi3(1) =", Phi3.subs(w,1))
print("D3 =", D3, " degree:", sp.Poly(D3, s, r).total_degree())
C3 = cusps_of(p3, Phi3)
print("cusps:", len(C3))
for wv, sv, rv in C3:
    real = abs(wv.imag) < 1e-30
    print(f"  t={wv:.5f} -> ({sv:.5f},{rv:.5f}) {'REAL' if real else 'cplx'}")
N3, gaps3, elim3 = nodes_of(p3, Phi3)
print(f"nodes: {len(N3)} (predict 1); intra-pair image gaps: {gaps3}")
for p_, q_ in N3:
    print(f"  pair t1={p_[0]:.5f} t2={q_[0]:.5f} -> ({p_[1]:.6f},{p_[2]:.6f})")
print(f"budget: {len(C3)} + {len(N3)} = {len(C3)+len(N3)}  (target 3 = (4-1)(4-2)/2)")
