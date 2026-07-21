"""
Note 7, stage A2: bitangent nodes of the F5 wall (the (2,2,1,1) fibers).
eq1 = (p(w2)-p(w1))/(w2-w1)  [deg 4]
eq2 = (Phi(w2)-Phi(w1))/(w2-w1) - p(w1)  [deg 5]
Expect 6 off-diagonal unordered pairs (genus bookkeeping), no triple points.
"""
import sympy as sp, json
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

w, w1, w2 = symbols("w w1 w2")
p5 = -w**5 + w**4 - R(14,5)*w**2 + R(9,5)*w
Phi5 = expand(integrate(p5, w))

eq1 = sp.expand(cancel((p5.subs(w, w2) - p5.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi5.subs(w, w2) - Phi5.subs(w, w1))/(w2 - w1) - p5.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(), flush=True)
print("lex GB over Q[w1,w2] ...", flush=True)
gb = groebner([eq1, eq2], w1, w2, order="lex")
print("GB sizes:", [sp.Poly(g.expr, w1, w2).total_degree() for g in gb.polys], flush=True)

elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
print("eliminant degree:", sp.Poly(elim, w2).degree())
with open("atlas5_bitangent_eliminant.txt", "w") as f:
    f.write(str(sp.factor(elim)))
print("eliminant factored:", sp.factor(elim), flush=True)

cands = sp.nroots(sp.Poly(elim, w2), n=50, maxsteps=800)
rel = [g.expr for g in gb.polys if set(g.free_symbols) == {w1, w2}]
print("linking polys:", len(rel))

def rnd(v, nd=9):
    c = complex(sp.N(v, 60))
    return complex(round(c.real, nd), round(c.imag, nd))

pairs = []
for cv in cands:
    rsets = []
    for g in rel:
        sub = sp.Poly(sp.expand(g.subs(w2, cv)), w1)
        rsets.append(set(rnd(v) for v in sp.nroots(sub, n=60, maxsteps=800)))
    if not rsets:
        continue
    common = rsets[0]
    for rr in rsets[1:]:
        common &= rr
    for wv1 in common:
        a, b = complex(wv1), complex(cv)
        e1v = abs(complex(eq1.subs({w1: a, w2: b})))
        e2v = abs(complex(eq2.subs({w1: a, w2: b})))
        if e1v < 1e-8 and e2v < 1e-8:
            pairs.append((a, b))
print(f"residual-filtered ordered pairs: {len(pairs)}", flush=True)

# de-duplicate unordered pairs; classify diagonal (cusp evaporation) vs true nodes
cusp_ts = [complex(v) for v in sp.nroots(diff(p5, w), n=50)]
nodes = []
seen = []
for a, b in pairs:
    if abs(a - b) < 1e-7:
        continue
    key = frozenset((round(a.real,6), round(a.imag,6), round(b.real,6), round(b.imag,6)))
    if key in seen:
        continue
    seen.append(key)
    sv = complex(p5.subs(w, a)); rv = complex(a*p5.subs(w, a) - Phi5.subs(w, a))
    nodes.append((a, b, sv, rv))
print(f"\nunordered off-diagonal bitangent pairs: {len(nodes)}")
for a, b, sv, rv in nodes:
    realline = abs(sv.imag) < 1e-12 and abs(rv.imag) < 1e-12
    realc = abs(a.imag) < 1e-12 and abs(b.imag) < 1e-12
    print(f"  t1={a:.6f}  t2={b:.6f} -> (s,r)=({sv:.8f},{rv:.8f})  line-real={realline} contacts-real={realc}", flush=True)

# triple points: same (s,r) hit by >=2 distinct pairs -> (2,2,2)
trip = []
for i in range(len(nodes)):
    for j in range(i+1, len(nodes)):
        if abs(nodes[i][2]-nodes[j][2]) < 1e-9 and abs(nodes[i][3]-nodes[j][3]) < 1e-9:
            trip.append((i, j))
print("\ntriple points ((2,2,2) fibers):", len(trip), flush=True)

# node == cusp point? -> (3,2,*,*) overlap
cusp_pts = []
for ct in cusp_ts:
    cusp_pts.append((complex(p5.subs(w, ct)), complex(ct*p5.subs(w, ct) - Phi5.subs(w, ct))))
overl = 0
for a, b, sv, rv in nodes:
    for csv, crv in cusp_pts:
        if abs(sv-csv) < 1e-9 and abs(rv-crv) < 1e-9:
            overl += 1
print("node-cusp overlaps ((3,2) strata):", overl)

# genus bookkeeping
delta = len(nodes)*1 + 4*1 + len(trip)*2 + overl*2     # triple pt costs 3 pairs => +2 extra, tacnode +2
print(f"\ngenus bookkeeping: 6 nodes?? {len(nodes)} + 4 cusps + extras = {delta}  (target (6-1)(6-2)/2 = 10)")

json.dump({"n_nodes": len(nodes),
           "nodes": [[str(a), str(b), repr(sv), repr(rv)] for a,b,sv,rv in nodes],
           "triple_points": len(trip), "node_cusp_overlaps": overl, "delta_sum": delta},
          open("atlas5_bitangents.json","w"), indent=1)
