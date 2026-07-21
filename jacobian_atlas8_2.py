"""
Note 10, stage B: bitangents/nodes of the F8 wall. Predict 21 nodes (budget 7+21=28).
Eliminant should factor: (12*p8')^2 * cofactor[deg 42], total deg 56 = (n-1)(n-2).
Protocol: split factors FIRST (lesson: never root across double roots), then cluster
the 42 contact images (s,r) at 110 digits into 21 pairs (note-9 scheme, mpmath-native).
Dumps: atlas8_bitangent_eliminant.txt, atlas8_elim_raw.txt, atlas8_cofactor.txt, atlas8_bitangents.json
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

t0 = time.time()
w, w1, w2 = symbols("w w1 w2")
p8 = -w**8 + w**7 - R(35,12)*w**2 + R(23,12)*w
Phi9 = expand(integrate(p8, w))

eq1 = sp.expand(cancel((p8.subs(w, w2) - p8.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi9.subs(w, w2) - Phi9.subs(w, w1))/(w2 - w1) - p8.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(),
      "  (Bezout 56 ordered)", flush=True)

gb = groebner([eq1, eq2], w1, w2, order="lex")
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
ec = sp.gcd_list([c for c in sp.Poly(elim, w2).coeffs()])
elim = sp.expand(elim/ec)
print(f"eliminant degree: {sp.Poly(elim, w2).degree()}  [predict 56]  [{time.time()-t0:.0f}s]", flush=True)
with open("atlas8_elim_raw.txt","w") as f: f.write(str(sp.expand(elim)))

p8p = diff(p8, w).subs(w, w2)
q1, rmd = sp.div(elim, p8p, w2); assert expand(rmd) == 0, "p8' does not divide eliminant!"
q2, rmd = sp.div(q1, p8p, w2); assert expand(rmd) == 0, "p8' does not divide TWICE!"
cofactor = sp.expand(q2 / 144)
scale_check = sp.expand(elim - 144*cofactor*p8p**2) == 0
print(f"elim = (12*p8')^2 * cofactor[deg {sp.degree(cofactor, w2)}] EXACT: {scale_check}"
      f"  (cofactor deg predict 42)", flush=True)
with open("atlas8_bitangent_eliminant.txt","w") as f:
    f.write("(12*p8')^2 * cofactor[deg42]; elim = 144*p8'^2*cofactor exactly\n")
with open("atlas8_cofactor.txt","w") as f: f.write(str(cofactor))

# --- exact transversality certificates over Q ---
g1 = sp.gcd(cofactor, diff(cofactor, w2))
g2 = sp.gcd(cofactor, p8p)
sqfree = (sp.expand(g1) in (sp.Integer(1), sp.Integer(-1)))
coprime = (sp.expand(g2) in (sp.Integer(1), sp.Integer(-1)))
print(f"[{time.time()-t0:.0f}s] cofactor squarefree over Q: {sqfree} | coprime to p8': {coprime}", flush=True)

# --- roots: factors separately (lesson), 110 digits, mpmath-native clustering ---
import mpmath as mp
mp.mp.dps = 110
print(f"[{time.time()-t0:.0f}s] rooting factors at 110 digits ...", flush=True)
def mpify(v): return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
cusp_roots = [mpify(v) for v in sp.nroots(diff(p8, w), n=110, maxsteps=4000)]
cof_roots  = [mpify(v) for v in sp.nroots(sp.Poly(cofactor, w2), n=110, maxsteps=10000)]
print(f"  cofactor roots: {len(cof_roots)}  [predict 42]", flush=True)

def mp_wall_pt(t):
    sv = -t**8 + t**7 - mp.mpf(35)/12*t**2 + mp.mpf(23)/12*t
    phit = -t**9/9 + t**8/8 - mp.mpf(35)/36*t**3 + mp.mpf(23)/24*t**2
    return sv, t*sv - phit
maxres = mp.mpf(0)
cofm = [mp.mpc(mp.mpf(str(sp.re(c))), mp.mpf(str(sp.im(c)))) for c in sp.Poly(cofactor, w2).all_coeffs()]
for v in cof_roots:
    val = mp.polyval(cofm, v)
    maxres = max(maxres, abs(val))
print(f"  max |cofactor(root)| = {mp.nstr(maxres, 4)}", flush=True)

imgs = [(v, *mp_wall_pt(v)) for v in cof_roots]
used = [False]*len(imgs); pairs = []
for i, (a, sa, ra) in enumerate(imgs):
    if used[i]: continue
    bd, bd2, best = mp.mpf("1e60"), mp.mpf("1e60"), -1
    for j in range(i+1, len(imgs)):
        if used[j]: continue
        ddst = abs(imgs[j][1]-sa) + abs(imgs[j][2]-ra)
        if ddst < bd: bd2, bd, best = bd, ddst, j
        elif ddst < bd2: bd2 = ddst
    assert best >= 0
    pairs.append((a, imgs[best][0], (sa+imgs[best][1])/2, (ra+imgs[best][2])/2, bd, bd2))
    used[i] = used[best] = True
worst_gap = max(p_[4] for p_ in pairs); best_margin = min(p_[5]/max(p_[4], mp.mpf("1e-999")) for p_ in pairs)
print(f"unordered bitangent pairs: {len(pairs)}  [predict 21]")
print(f"  worst intra-pair gap: {mp.nstr(worst_gap,4)}  |  min margin (2nd-best/best): {mp.nstr(best_margin,4)}  (must be huge)")

nodes = []
def res_eq(eqq, a, b):
    v = eqq.subs({w1: complex(a), w2: complex(b)})
    return abs(complex(sp.N(v, 30)))
for a, b, sv, rv, gap, gap2 in pairs:
    mx = max(res_eq(eq1, a, b), res_eq(eq2, a, b), res_eq(eq1, b, a), res_eq(eq2, b, a))
    realline = abs(mp.im(sv)) < mp.mpf("1e-40") and abs(mp.im(rv)) < mp.mpf("1e-40")
    realc = abs(mp.im(a)) < mp.mpf("1e-60") and abs(mp.im(b)) < mp.mpf("1e-60")
    conj = abs(mp.conj(a) - b) < mp.mpf("1e-60")
    kind = "CRUNODE" if (realline and realc) else ("ACNODE" if realline else "complex")
    nodes.append((a, b, sv, rv, kind, mx, gap))
print(f"  worst eq-residual: {max(n_[5] for n_ in nodes):.3e}")
n_cr = sum(1 for n_ in nodes if n_[4] == "CRUNODE"); n_ac = sum(1 for n_ in nodes if n_[4] == "ACNODE")
print(f"CRUNODES: {n_cr} [predict 0]   ACNODES: {n_ac} [guess 3]   complex: {len(nodes)-n_cr-n_ac}")
for a, b, sv, rv, kind, mx, gap in nodes:
    if kind != "complex":
        print(f"  t1={mp.nstr(a,10)} t2={mp.nstr(b,10)} -> ({mp.nstr(sv,12)}, {mp.nstr(rv,12)})  {kind}", flush=True)

trip = sum(1 for i in range(len(nodes)) for j in range(i+1, len(nodes))
           if abs(nodes[i][2]-nodes[j][2]) < mp.mpf("1e-40") and abs(nodes[i][3]-nodes[j][3]) < mp.mpf("1e-40"))
cusp_pts = [mp_wall_pt(c) for c in cusp_roots]
overl = sum(1 for a,b,sv,rv,k,m,g in nodes for csv,crv in cusp_pts
            if abs(sv-csv) < mp.mpf("1e-40") and abs(rv-crv) < mp.mpf("1e-40"))
print(f"triple points: {trip}   node-cusp overlaps: {overl}")
delta = len(nodes) + 7 + trip*2 + overl*2
print(f"delta budget: {len(nodes)} nodes + 7 cusps = {delta}  (target 28)  balanced: {delta==28}")
json.dump({"n_nodes": len(nodes), "crunodes": n_cr, "acnodes": n_ac, "delta": delta,
           "eliminant_scale": bool(scale_check), "cofactor_squarefree": bool(sqfree),
           "cofactor_coprime_p8p": bool(coprime), "worst_pair_gap": str(worst_gap),
           "min_margin": str(best_margin),
           "nodes": [[mp.nstr(a,30), mp.nstr(b,30), mp.nstr(sv,30), mp.nstr(rv,30), k] for a,b,sv,rv,k,m,g in nodes]},
          open("atlas8_bitangents.json","w"), indent=1, default=str)
print(f"[stage B done {time.time()-t0:.0f}s]", flush=True)
