"""
Note 11, stage B: 28 nodes of the decic wall (budget 8 + 28 = 36).
Predict: eliminant = (15*p9')^2 * cofactor[deg 56]; cofactor 56 roots -> 28 pairs;
1 crunode (dance 0,1,0,1,0 -> 1), 3 acnodes (staircase 1,1,2,2,3 -> 3).
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

t0 = time.time()
w, w1, w2 = symbols("w w1 w2")
p9 = -w**9 + w**8 - R(44,15)*w**2 + R(29,15)*w
Phi10 = expand(integrate(p9, w))

eq1 = sp.expand(cancel((p9.subs(w, w2) - p9.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi10.subs(w, w2) - Phi10.subs(w, w1))/(w2 - w1) - p9.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(),
      "  (Bezout 72 ordered)", flush=True)

gb = groebner([eq1, eq2], w1, w2, order="lex")
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
ec = sp.gcd_list([c for c in sp.Poly(elim, w2).coeffs()])
elim = sp.expand(elim/ec)
print(f"eliminant degree: {sp.Poly(elim, w2).degree()}  [predict 72]  [{time.time()-t0:.0f}s]", flush=True)
with open("atlas9_elim_raw.txt","w") as f: f.write(str(elim))

p9p = diff(p9, w).subs(w, w2)
q1, rmd = sp.div(elim, p9p, w2); assert expand(rmd) == 0
q2, rmd = sp.div(q1, p9p, w2); assert expand(rmd) == 0
cofactor = sp.expand(q2 / 225)          # 15^2 = 225
scale_check = sp.expand(elim - 225*cofactor*p9p**2) == 0
print(f"elim = (15*p9')^2 * cofactor[deg {sp.degree(cofactor, w2)}] EXACT: {scale_check} (predict 56)", flush=True)
with open("atlas9_cofactor.txt","w") as f: f.write(str(cofactor))

g1 = sp.gcd(cofactor, diff(cofactor, w2)); g2 = sp.gcd(cofactor, p9p)
sqfree = (sp.degree(g1, w2) == 0); coprime = (sp.degree(g2, w2) == 0)
print(f"[{time.time()-t0:.0f}s] cofactor squarefree over Q: {sqfree} | coprime to p9': {coprime}", flush=True)

import mpmath as mp
mp.mp.dps = 110
def mpify(v): return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
cusp_roots = [mpify(v) for v in sp.nroots(diff(p9, w), n=110, maxsteps=8000)]
cof_roots  = [mpify(v) for v in sp.nroots(sp.Poly(cofactor, w2), n=110, maxsteps=15000)]
print(f"  cofactor roots: {len(cof_roots)}  [predict 56]", flush=True)

def mp_wall_pt(t):
    sv = -t**9 + t**8 - mp.mpf(44)/15*t**2 + mp.mpf(29)/15*t
    phit = -t**10/10 + t**9/9 - mp.mpf(44)/45*t**3 + mp.mpf(29)/30*t**2
    return sv, t*sv - phit
cofm = [mp.mpc(mp.mpf(str(sp.re(c))), mp.mpf(str(sp.im(c)))) for c in sp.Poly(cofactor, w2).all_coeffs()]
maxres = max(abs(mp.polyval(cofm, v)) for v in cof_roots)
print(f"  max |cofactor(root)| = {mp.nstr(maxres, 4)}", flush=True)

imgs = [(v, *mp_wall_pt(v)) for v in cof_roots]
used = [False]*len(imgs); pairs = []
for i, (a, sa, ra) in enumerate(imgs):
    if used[i]: continue
    bd, bd2, best = mp.mpf("1e60"), mp.mpf("1e60"), -1
    for j in range(i+1, len(imgs)):
        if used[j]: continue
        dd = abs(imgs[j][1]-sa) + abs(imgs[j][2]-ra)
        if dd < bd: bd2, bd, best = bd, dd, j
        elif dd < bd2: bd2 = dd
    assert best >= 0
    pairs.append((a, imgs[best][0], (sa+imgs[best][1])/2, (ra+imgs[best][2])/2, bd, bd2))
    used[i] = used[best] = True
worst_gap = max(p_[4] for p_ in pairs); best_margin = min(p_[5]/max(p_[4], mp.mpf("1e-999")) for p_ in pairs)
print(f"unordered bitangent pairs: {len(pairs)}  [predict 28]")
print(f"  worst intra-pair gap: {mp.nstr(worst_gap,4)}  |  min margin: {mp.nstr(best_margin,4)}")

def res_eq(eqq, a, b):
    return abs(complex(sp.N(eqq.subs({w1: complex(a), w2: complex(b)}), 30)))
nodes = []
for a, b, sv, rv, gap, gap2 in pairs:
    mx = max(res_eq(eq1, a, b), res_eq(eq2, a, b), res_eq(eq1, b, a), res_eq(eq2, b, a))
    realline = abs(mp.im(sv)) < mp.mpf("1e-40") and abs(mp.im(rv)) < mp.mpf("1e-40")
    realc = abs(mp.im(a)) < mp.mpf("1e-60") and abs(mp.im(b)) < mp.mpf("1e-60")
    kind = "CRUNODE" if (realline and realc) else ("ACNODE" if realline else "complex")
    nodes.append((a, b, sv, rv, kind, mx, gap))
print(f"  worst eq-residual: {max(n_[5] for n_ in nodes):.3e}")
n_cr = sum(1 for n_ in nodes if n_[4] == "CRUNODE"); n_ac = sum(1 for n_ in nodes if n_[4] == "ACNODE")
print(f"CRUNODES: {n_cr} [predict 1]   ACNODES: {n_ac} [predict 3]   complex: {len(nodes)-n_cr-n_ac}")
for a, b, sv, rv, kind, mx, gap in nodes:
    if kind != "complex":
        print(f"  t1={mp.nstr(a,10)} t2={mp.nstr(b,10)} -> ({mp.nstr(sv,12)}, {mp.nstr(rv,12)})  {kind}", flush=True)

# crunode whisker-modality: fiber at crunode - real residual roots?
for a, b, sv, rv, kind, mx, gap in nodes:
    if kind == "CRUNODE":
        import numpy as np
        hco = [-1/10, 1/9, 0,0,0,0,0, -44/45, 29/30, -float(mp.re(sv)), float(mp.re(rv))]
        roots = np.roots([complex(x) for x in hco])
        others = [zz for zz in roots if abs(zz - complex(a)) > 1e-4 and abs(zz - complex(b)) > 1e-4]
        nre = sum(1 for zz in others if abs(zz.imag) < 1e-8)
        nre_esc = sum(1 for zz in others if abs(zz.imag) < 1e-8 and abs(float(mp.re(sv)) - (-(zz.real**9)+zz.real**8-44*zz.real**2/15+29*zz.real/15)) < 1e-6)
        print(f"  crunode modality: {nre} real residual roots, of which {nre_esc} with gamma~0 (escaping)")

trip = sum(1 for i in range(len(nodes)) for j in range(i+1, len(nodes))
           if abs(nodes[i][2]-nodes[j][2]) < mp.mpf("1e-40") and abs(nodes[i][3]-nodes[j][3]) < mp.mpf("1e-40"))
cusp_pts = [mp_wall_pt(c) for c in cusp_roots]
overl = sum(1 for a,b,sv,rv,k,m,g in nodes for csv,crv in cusp_pts
            if abs(sv-csv) < mp.mpf("1e-40") and abs(rv-crv) < mp.mpf("1e-40"))
print(f"triple points: {trip}   node-cusp overlaps: {overl}")
delta = len(nodes) + 8 + trip*2 + overl*2
print(f"delta budget: {len(nodes)} nodes + 8 cusps = {delta}  (target 36)  balanced: {delta==36}")
json.dump({"n_nodes": len(nodes), "crunodes": n_cr, "acnodes": n_ac, "delta": delta,
           "eliminant_scale": bool(scale_check), "cofactor_squarefree": bool(sqfree),
           "cofactor_coprime_p9p": bool(coprime), "worst_pair_gap": str(worst_gap),
           "nodes": [[mp.nstr(a,30), mp.nstr(b,30), mp.nstr(sv,30), mp.nstr(rv,30), k] for a,b,sv,rv,k,m,g in nodes]},
          open("atlas9_bitangents.json","w"), indent=1, default=str)
print(f"[stage B done {time.time()-t0:.0f}s]", flush=True)
