import sympy as sp
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

w, w1, w2 = symbols("w w1 w2")
p4 = -w**4 + w**3 - R(27,10)*w**2 + R(17,10)*w
Phi4 = expand(integrate(p4, w))

eq1 = sp.expand(cancel((p4.subs(w, w2) - p4.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi4.subs(w, w2) - Phi4.subs(w, w1))/(w2 - w1) - p4.subs(w, w1)))
gb = groebner([eq1, eq2], w1, w2, order="lex")
print("GB polys:", [sp.Poly(g.expr, w1, w2).total_degree() for g in gb.polys])
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
print("eliminant:", sp.factor(elim))
cands = sp.nroots(sp.Poly(elim, w2), n=40, maxsteps=500)
rel = [g.expr for g in gb.polys if set(g.free_symbols) == {w1, w2}]
print("linking poly degrees:", [sp.Poly(g, w1, w2).total_degree() for g in rel])

pairs = []
for cv in cands:
    # solve linking polys for w1 at w2 = cv: intersect root sets at high precision
    rsets = []
    for g in rel:
        sub = sp.Poly(sp.expand(g.subs(w2, cv)), w1)
        def rnd(v):
            c = complex(sp.N(v, 40))
            return complex(round(c.real, 8), round(c.imag, 8))
        rsets.append(set(rnd(v) for v in sp.nroots(sub, n=40, maxsteps=500)))
    # also require original eqs to vanish
    common = rsets[0]
    for rr in rsets[1:]:
        common &= rr
    for wv1 in common:
        a, b = complex(wv1), complex(cv)
        # exact residual check on eq1, eq2 at (a,b)
        e1v = abs(complex(eq1.subs({w1: a, w2: b})))
        e2v = abs(complex(eq2.subs({w1: a, w2: b})))
        if e1v < 1e-10 and e2v < 1e-10:
            pairs.append((a, b))

print(f"residual-filtered pairs: {len(pairs)}")
for a, b in pairs:
    sv = complex(p4.subs(w, a)); rv = complex(a*p4.subs(w, a) - Phi4.subs(w, a))
    diag = abs(a-b) < 1e-7
    print(f"  w1={a:.7f} w2={b:.7f} {'[DIAG/cusp]' if diag else ''}")
    if not diag:
        realC = abs(a.imag) < 1e-10 and abs(b.imag) < 1e-10
        print(f"     -> (s,r) = ({sv:.8f}, {rv:.8f})  real contacts: {realC}, real line: {abs(sv.imag)<1e-10 and abs(rv.imag)<1e-10}")
