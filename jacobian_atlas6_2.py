"""
Note 8, stage A2: bitangents of the F6 wall. Predict 10 nodes (budget 5 + 10 = 15).
Contact eliminant should factor as p6'^2 * (degree-20 bitangent eliminant).
Also: triple points ((2,2,2)), node-cusp overlaps ((3,2)), support-diff of D4/D5/D6.
"""
import sympy as sp, json
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

w, w1, w2 = symbols("w w1 w2")
p6 = -w**6 + w**5 - R(20,7)*w**2 + R(13,7)*w
Phi6 = expand(integrate(p6, w))

eq1 = sp.expand(cancel((p6.subs(w, w2) - p6.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi6.subs(w, w2) - Phi6.subs(w, w1))/(w2 - w1) - p6.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(), flush=True)
gb = groebner([eq1, eq2], w1, w2, order="lex")
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
print("eliminant degree:", sp.Poly(elim, w2).degree(), flush=True)
fel = sp.factor(elim)
print("eliminant factored:", fel, flush=True)
with open("atlas6_bitangent_eliminant.txt","w") as f: f.write(str(fel))

cands = sp.nroots(sp.Poly(elim, w2), n=50, maxsteps=1500)
rel = [g.expr for g in gb.polys if set(g.free_symbols) == {w1, w2}]

def rnd(v, nd=9):
    c = complex(sp.N(v, 60))
    return complex(round(c.real, nd), round(c.imag, nd))

pairs = []
for cv in cands:
    rsets = []
    for g in rel:
        sub = sp.Poly(sp.expand(g.subs(w2, cv)), w1)
        rsets.append(set(rnd(v) for v in sp.nroots(sub, n=60, maxsteps=1000)))
    if not rsets:
        continue
    common = rsets[0]
    for rr in rsets[1:]:
        common &= rr
    for wv1 in common:
        a, b = complex(wv1), complex(cv)
        if abs(complex(eq1.subs({w1: a, w2: b}))) < 1e-8 and abs(complex(eq2.subs({w1: a, w2: b}))) < 1e-8:
            pairs.append((a, b))
print(f"residual-filtered ordered solutions: {len(pairs)} (Bezout 30)", flush=True)

nodes = []
seen = []
for a, b in pairs:
    if abs(a - b) < 1e-7:
        continue
    key = frozenset((round(a.real,6), round(a.imag,6), round(b.real,6), round(b.imag,6)))
    if key in seen:
        continue
    seen.append(key)
    sv = complex(p6.subs(w, a)); rv = complex(a*p6.subs(w, a) - Phi6.subs(w, a))
    nodes.append((a, b, sv, rv))
print(f"unordered off-diagonal bitangent pairs: {len(nodes)}  [predict 10]", flush=True)
for a, b, sv, rv in nodes:
    realline = abs(sv.imag) < 1e-10 and abs(rv.imag) < 1e-10
    realc = abs(a.imag) < 1e-10 and abs(b.imag) < 1e-10
    kind = "CRUNODE" if realline and realc else ("ACNODE" if realline else "complex")
    print(f"  t1={a:.6f} t2={b:.6f} -> ({sv:.8f},{rv:.8f})  {kind}", flush=True)

trip = sum(1 for i in range(len(nodes)) for j in range(i+1, len(nodes))
           if abs(nodes[i][2]-nodes[j][2]) < 1e-9 and abs(nodes[i][3]-nodes[j][3]) < 1e-9)
print("triple points ((2,2,2)):", trip, flush=True)

cusp_ts = [complex(v) for v in sp.nroots(diff(p6, w), n=50)]
cusp_pts = [(complex(p6.subs(w, ct)), complex(ct*p6.subs(w, ct) - Phi6.subs(w, ct))) for ct in cusp_ts]
overl = sum(1 for a,b,sv,rv in nodes for csv,crv in cusp_pts
            if abs(sv-csv) < 1e-9 and abs(rv-crv) < 1e-9)
print("node-cusp overlaps ((3,2) strata):", overl, flush=True)

delta = len(nodes) + 5 + trip*2 + overl*2
print(f"delta budget: {len(nodes)} nodes + 5 cusps = {delta}  (target 15)  balanced: {delta==15}")

json.dump({"n_nodes": len(nodes),
           "nodes": [[str(a), str(b), repr(sv), repr(rv)] for a,b,sv,rv in nodes],
           "triple_points": trip, "node_cusp_overlaps": overl, "delta_sum": delta},
          open("atlas6_bitangents.json","w"), indent=1)

# support diff D4/D5/D6
def supp(path_terms):
    return set((e1, e2) for _, e1, e2 in path_terms)
D4 = [(-20000000,4,0),(20000000,3,1),(18900000,3,0),(-57100000,2,2),(75376000,2,1),
      (-46613461,2,0),(38210000,1,3),(-67741050,1,2),(45893931,1,1),(-6657115,1,0),
      (-8192000,0,5),(17058675,0,4),(-12278715,0,3),(1957975,0,2)]
D5 = [(3037500000,5,0),(-3037500000,4,1),(11907000000,4,0),(-50625000,3,2),
      (-19265850000,3,1),(15615855000,3,0),(9797625000,2,3),(8794980000,2,2),
      (-23436459000,2,1),(8496467856,2,0),(-6968981250,1,4),(2780676000,1,3),
      (8071107300,1,2),(-6147920736,1,1),(892142910,1,0),(1220703125,0,6),
      (-46159500,0,5),(-2098753050,0,4),(1564980260,0,3),(-247817475,0,2)]
D6p = sp.Poly(sp.sympify(open("atlas6_wall.txt").read()), sp.symbols("s r"))
s, r = sp.symbols("s r")
D6s = set((exp[s], exp[r]) for exp in sp.Poly(D6p.as_expr(), s, r).monoms())
S4, S5 = supp(D4), supp(D5)
for name, S_, n in [("D4", S4, 5), ("D5", S5, 6), ("D6", D6s, 7)]:
    full = set((i, j) for i in range(n) for j in range(n+1) if 2 <= i+j <= n) | set()
    allmon = set((i, j) for i in range(0, n) for j in range(0, n+1) if i+j <= n and (i, j) != (0,0))
    missing = sorted(allmon - S_, key=lambda e: (e[0]+e[1], e))
    print(f"{name}: {len(S_)} terms; missing from full triangle deg<= {n}: {missing}")
