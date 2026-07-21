"""
Note 9, stage A2: bitangents of the F7 wall. Predict 15 nodes (budget 6 + 15 = 21).
Contact eliminant should be (c*p7')^2 * (degree-30 bitangent eliminant), Bezout 42.
"""
import sympy as sp, json
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

w, w1, w2 = symbols("w w1 w2")
p7 = -w**7 + w**6 - R(81,28)*w**2 + R(53,28)*w
Phi7 = expand(integrate(p7, w))

eq1 = sp.expand(cancel((p7.subs(w, w2) - p7.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi7.subs(w, w2) - Phi7.subs(w, w1))/(w2 - w1) - p7.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(), flush=True)
gb = groebner([eq1, eq2], w1, w2, order="lex")
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
print("eliminant degree:", sp.Poly(elim, w2).degree(), flush=True)
fel = sp.factor(elim)
print("eliminant factored:", str(fel)[:600], "...", flush=True)
with open("atlas7_bitangent_eliminant.txt","w") as f: f.write(str(fel))
# check squared factor = multiple of p7'
p7p = diff(p7, w)
if fel.is_Mul:
    for fac in fel.args:
        if fac.is_Pow and fac.exp == 2 or (fac.is_Poly):
            pass
sq = None
for fac in (fel.args if fel.is_Mul else [fel]):
    base, exp_ = (fac.args if fac.is_Pow else (fac, 1))
    if exp_ == 2 and sp.degree(base, w2) == 6:
        sq = base
if sq is not None:
    ratio = sp.cancel(sq / p7p.subs(w, w2))
    print("squared sextic factor / p7' =", ratio, flush=True)

cands = sp.nroots(sp.Poly(elim, w2), n=50, maxsteps=2500)
rel = [g.expr for g in gb.polys if set(g.free_symbols) == {w1, w2}]

def rnd(v, nd=9):
    c = complex(sp.N(v, 60))
    return complex(round(c.real, nd), round(c.imag, nd))

pairs = []
for cv in cands:
    rsets = []
    for g in rel:
        sub = sp.Poly(sp.expand(g.subs(w2, cv)), w1)
        rsets.append(set(rnd(v) for v in sp.nroots(sub, n=60, maxsteps=1600)))
    if not rsets:
        continue
    common = rsets[0]
    for rr in rsets[1:]:
        common &= rr
    for wv1 in common:
        a, b = complex(wv1), complex(cv)
        if abs(complex(eq1.subs({w1: a, w2: b}))) < 1e-8 and abs(complex(eq2.subs({w1: a, w2: b}))) < 1e-8:
            pairs.append((a, b))
print(f"residual-filtered ordered solutions: {len(pairs)} (Bezout 42)", flush=True)

nodes = []
seen = []
for a, b in pairs:
    if abs(a - b) < 1e-7:
        continue
    key = frozenset((round(a.real,6), round(a.imag,6), round(b.real,6), round(b.imag,6)))
    if key in seen:
        continue
    seen.append(key)
    sv = complex(p7.subs(w, a)); rv = complex(a*p7.subs(w, a) - Phi7.subs(w, a))
    nodes.append((a, b, sv, rv))
print(f"unordered bitangent pairs: {len(nodes)}  [predict 15]", flush=True)
n_crunode = n_acnode = 0
for a, b, sv, rv in nodes:
    realline = abs(sv.imag) < 1e-9 and abs(rv.imag) < 1e-9
    realc = abs(a.imag) < 1e-10 and abs(b.imag) < 1e-10
    kind = "CRUNODE" if realline and realc else ("ACNODE" if realline else "complex")
    if kind == "CRUNODE": n_crunode += 1
    if kind == "ACNODE": n_acnode += 1
    if realline:
        print(f"  t1={a:.6f} t2={b:.6f} -> ({sv:.8f},{rv:.8f})  {kind}", flush=True)
print(f"real nodes: {n_crunode} crunodes + {n_acnode} acnodes", flush=True)

trip = sum(1 for i in range(len(nodes)) for j in range(i+1, len(nodes))
           if abs(nodes[i][2]-nodes[j][2]) < 1e-9 and abs(nodes[i][3]-nodes[j][3]) < 1e-9)
print("triple points:", trip, flush=True)
cusp_ts = [complex(v) for v in sp.nroots(diff(p7, w), n=50)]
cusp_pts = [(complex(p7.subs(w, ct)), complex(ct*p7.subs(w, ct) - Phi7.subs(w, ct))) for ct in cusp_ts]
overl = sum(1 for a,b,sv,rv in nodes for csv,crv in cusp_pts
            if abs(sv-csv) < 1e-9 and abs(rv-crv) < 1e-9)
print("node-cusp overlaps:", overl, flush=True)

delta = len(nodes) + 6 + trip*2 + overl*2
print(f"delta budget: {len(nodes)} + 6 = {delta}  (target 21)  balanced: {delta==21}")
json.dump({"n_nodes": len(nodes), "crunodes": n_crunode, "acnodes": n_acnode,
           "nodes": [[str(a), str(b), repr(sv), repr(rv)] for a,b,sv,rv in nodes],
           "triple_points": trip, "node_cusp_overlaps": overl, "delta_sum": delta},
          open("atlas7_bitangents.json","w"), indent=1)
