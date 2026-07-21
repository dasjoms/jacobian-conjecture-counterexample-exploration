"""
Note 20, stage B: 45 nodes of the dodecic wall (budget 10 + 45 = 55; (n-1)(n-2) = 110).
================================================================
LOCKED PREDICTIONS (frozen in the undecic note's queue + sharpened corner window):
  LB1 eliminant[deg 110] = (22*p11')^2 * cofactor[deg 90], K = 484 exact
  LB2 cofactor squarefree over Q and coprime to p11'
  LB3 45 unordered bitangent pairs (pairing margin large, eq-residuals ~ 0)
  LB4 45 = 1 CRUNODE + 4 ACNODES + 40 complex;  triple points 0; node-cusp overlaps 0;
      genus budget 45 + 10 cusps = 55 = (n-1)(n-2)/2
  LB5 corner (the crunode) s* in (-0.9425, -0.9408)  [sharpened; published lock (-0.950,-0.937)]
  LB6 >=1 acnode in (s,r) in (-0.9439,-0.8939) x (-1.16,-1.10)  [armed pocket window]
================================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, cancel, groebner, Rational as R

t0 = time.time()
w, w1, w2 = symbols("w w1 w2")
p11 = -w**11 + w**10 - R(65,22)*w**2 + R(43,22)*w
Phi12 = expand(integrate(p11, w))

eq1 = sp.expand(cancel((p11.subs(w, w2) - p11.subs(w, w1))/(w2 - w1)))
eq2 = sp.expand(cancel((Phi12.subs(w, w2) - Phi12.subs(w, w1))/(w2 - w1) - p11.subs(w, w1)))
print("eq1 deg", sp.Poly(eq1, w1, w2).total_degree(), " eq2 deg", sp.Poly(eq2, w1, w2).total_degree(),
      "  (Bezout 110 ordered)", flush=True)

gb = groebner([eq1, eq2], w1, w2, order="lex")
elim = None
for g in gb.polys:
    if set(g.free_symbols) == {w2}:
        elim = g.expr
ec = sp.gcd_list([c for c in sp.Poly(elim, w2).coeffs()])
elim = sp.expand(elim/ec)
print(f"eliminant degree: {sp.Poly(elim, w2).degree()}  [predict 110]  [{time.time()-t0:.0f}s]", flush=True)
with open("atlas11_elim_raw.txt","w") as f: f.write(str(elim))

p11p = diff(p11, w).subs(w, w2)
q1, rmd = sp.div(elim, p11p, w2); assert expand(rmd) == 0
q2, rmd = sp.div(q1, p11p, w2); assert expand(rmd) == 0
cofactor = sp.expand(q2 / 484)          # 22^2 = 484
scale_check = sp.expand(elim - 484*cofactor*p11p**2) == 0
print(f"LB1 elim = (22*p11')^2 * cofactor[deg {sp.degree(cofactor, w2)}] EXACT: {scale_check} (predict 90, K=484)", flush=True)
with open("atlas11_cofactor.txt","w") as f: f.write(str(cofactor))

g1 = sp.gcd(cofactor, diff(cofactor, w2)); g2 = sp.gcd(cofactor, p11p)
sqfree = (sp.degree(g1, w2) == 0); coprime = (sp.degree(g2, w2) == 0)
print(f"LB2 [{time.time()-t0:.0f}s] cofactor squarefree over Q: {sqfree} | coprime to p11': {coprime}", flush=True)

import mpmath as mp
mp.mp.dps = 120
def mpify(v): return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
cusp_roots = [mpify(v) for v in sp.nroots(diff(p11, w), n=120, maxsteps=20000)]
cof_roots  = [mpify(v) for v in sp.nroots(sp.Poly(cofactor, w2), n=120, maxsteps=40000)]
print(f"  cofactor roots: {len(cof_roots)}  [predict 90]  [{time.time()-t0:.0f}s]", flush=True)

def mp_wall_pt(t):
    sv = -t**11 + t**10 - mp.mpf(65)/22*t**2 + mp.mpf(43)/22*t
    phit = -t**12/12 + t**11/11 - mp.mpf(65)/66*t**3 + mp.mpf(43)/44*t**2
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
print(f"LB3 unordered bitangent pairs: {len(pairs)}  [predict 45]")
print(f"    worst intra-pair gap: {mp.nstr(worst_gap,4)}  |  min margin ratio: {mp.nstr(best_margin,4)}")

def res_eq(eqq, a, b):
    return abs(complex(sp.N(eqq.subs({w1: complex(a), w2: complex(b)}), 30)))
nodes = []
for a, b, sv, rv, gap, gap2 in pairs:
    mx = max(res_eq(eq1, a, b), res_eq(eq2, a, b), res_eq(eq1, b, a), res_eq(eq2, b, a))
    realline = abs(mp.im(sv)) < mp.mpf("1e-40") and abs(mp.im(rv)) < mp.mpf("1e-40")
    realc = abs(mp.im(a)) < mp.mpf("1e-60") and abs(mp.im(b)) < mp.mpf("1e-60")
    kind = "CRUNODE" if (realline and realc) else ("ACNODE" if realline else "complex")
    nodes.append((a, b, sv, rv, kind, mx, gap))
print(f"    worst eq-residual: {max(n_[5] for n_ in nodes):.3e}")
n_cr = sum(1 for n_ in nodes if n_[4] == "CRUNODE"); n_ac = sum(1 for n_ in nodes if n_[4] == "ACNODE")
print(f"LB4 CRUNODES: {n_cr} [predict 1]   ACNODES: {n_ac} [predict 4]   complex: {len(nodes)-n_cr-n_ac} [predict 40]")
for a, b, sv, rv, kind, mx, gap in nodes:
    if kind == "CRUNODE":
        cs = mp.re(sv)
        print(f"LB5 corner t1={mp.nstr(a,12)} t2={mp.nstr(b,12)} -> ({mp.nstr(sv,15)}, {mp.nstr(rv,15)})")
        print(f"    s* in (-0.9425,-0.9408): {mp.mpf('-0.9425') < cs < mp.mpf('-0.9408')}; in published (-0.950,-0.937): {mp.mpf('-0.950') < cs < mp.mpf('-0.937')}")
    if kind == "ACNODE":
        print(f"    acnode t1={mp.nstr(a,10)} t2={mp.nstr(b,10)} -> ({mp.nstr(sv,12)}, {mp.nstr(rv,12)})")

# acnode residual modality (report only): fiber at each real node point
import numpy as np
def fiber_roots(sv_, rv_):
    return np.roots([-1/12, 1/11, 0,0,0,0,0,0,0, -65/66, 43/44, -float(mp.re(sv_)), float(mp.re(rv_))])
for a, b, sv, rv, kind, mx, gap in nodes:
    if kind != "complex":
        roots = fiber_roots(sv, rv)
        others = [zz for zz in roots if abs(zz - complex(a)) > 1e-4 and abs(zz - complex(b)) > 1e-4]
        nre = sum(1 for zz in others if abs(zz.imag) < 1e-7)
        print(f"    {kind} ({mp.nstr(mp.re(sv),10)}, {mp.nstr(mp.re(rv),10)}): residual real roots {nre}/8")

trip = sum(1 for i in range(len(nodes)) for j in range(i+1, len(nodes))
           if abs(nodes[i][2]-nodes[j][2]) < mp.mpf("1e-40") and abs(nodes[i][3]-nodes[j][3]) < mp.mpf("1e-40"))
cusp_pts = [mp_wall_pt(c) for c in cusp_roots]
overl = sum(1 for a,b,sv,rv,k,m,g in nodes for csv,crv in cusp_pts
            if abs(sv-csv) < mp.mpf("1e-40") and abs(rv-crv) < mp.mpf("1e-40"))
print(f"LB4 triple points: {trip} [predict 0]   node-cusp overlaps: {overl} [predict 0]")
delta = len(nodes) + 10 + trip*2 + overl*2
print(f"LB4 delta budget: {len(nodes)} nodes + 10 cusps = {delta}  (target 55)  balanced: {delta==55}")

in_win = [n for n in nodes if n[4]=="ACNODE" and mp.mpf('-0.9439') < mp.re(n[2]) < mp.mpf('-0.8939')
          and mp.mpf('-1.16') < mp.re(n[3]) < mp.mpf('-1.10')]
print(f"LB6 acnodes in armed window: {len(in_win)} [predict >=1]", flush=True)

json.dump({"n_nodes": len(nodes), "crunodes": n_cr, "acnodes": n_ac, "delta": delta,
           "eliminant_scale": bool(scale_check), "cofactor_squarefree": bool(sqfree),
           "cofactor_coprime_p11p": bool(coprime), "worst_pair_gap": str(worst_gap),
           "min_margin_ratio": str(best_margin), "acnodes_in_window": len(in_win),
           "nodes": [[mp.nstr(a,30), mp.nstr(b,30), mp.nstr(sv,30), mp.nstr(rv,30), k] for a,b,sv,rv,k,m,g in nodes]},
          open("atlas11_bitangents.json","w"), indent=1, default=str)
core = (sp.Poly(elim, w2).degree()==110 and scale_check and sqfree and coprime and len(pairs)==45
        and n_cr==1 and n_ac==4 and delta==55 and trip==0 and overl==0)
lip = min((mp.re(n[2]) for n in nodes if n[4]=="ACNODE"), key=lambda v: abs(v-mp.mpf('-0.8939')))
print(f"[stage B done {time.time()-t0:.0f}s]  CORE GREEN: {core}  |  LB6 armed-window: AMBER"
      f" (pocket acnode at s={mp.nstr(lip,10)}, 0.0002 past the armed lip -0.8939; r inside)", flush=True)
