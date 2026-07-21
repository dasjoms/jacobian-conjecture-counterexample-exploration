"""
Note 20, stage A: wall of F11 (fiber 12, d = 11, EVEN-fiber chamber - the CONE returns).
p11(w) = -w^11 + w^10 - (65/22)w^2 + (43/22)w ;  Phi12 = int p11
================================================================
LOCKED PREDICTIONS (frozen in the undecic note's queue + today's sharpenings):
  LA1 wall D11  : degree 12, irreducible, EXACTLY 76 TERMS = cone 79 - 3 genuine holes
                  {const, s^1, r^10}; (12,0) in supp (s^12 one weight step OUTSIDE the cone),
                  (0,12) fictitious; supp law n(n+1)/2+1-3 at n=12
  LA2           : D11(0,0) = 0 and s^2 | D11(s,0)          (from w^2 | Phi12)
  LA3 s^12 raw  : magnitude law 132^23 * 11^11 * (1/12)^11  (L = 132 = n(n-1)); ratio exactly 1
  LA4           : primitive [s^12] sign (+) = (-1)^12 (fingerprint 10/10);
                  content support subset of primes(12*11) = {2,3,11}
  LA5           : param identity t -> (p11(t), t*p11(t)-Phi12(t)) on wall, exact
  LA6 cusps     : gcd(p11',p11'') = 1 -> 10 simple ordinary cusps; Sturm: EXACTLY 2 real
                  (dance 1,2,1,2,...: odd d); BOTH real cusps HIT with residual nonic 9/9 real
  LA7 recipe    : kappa = -5 + 6/132 = -109/22; a = -(1+kappa)/(2+kappa) = -87/65;
                  den(p11') = 22, K = 484 = 22^2
================================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p11 = -w**11 + w**10 - R(65,22)*w**2 + R(43,22)*w
Phi12 = expand(integrate(p11, w))
assert expand(p11 - (2*w - 3*w**2 + w*(1-w)*(w**9 - R(1,22)))) == 0
assert p11.subs(w,1) == -1 and Phi12.subs(w,1) == 0
kappa = sp.diff(p11,w).subs(w,1)
print(f"p11   = {p11}\nPhi12 = {Phi12}")
print(f"kappa = {kappa} [lock -109/22] | recipe a = {-(1+kappa)/(2+kappa)} [lock -87/65]", flush=True)
out["kappa"] = str(kappa); out["a_recipe"] = str(-(1+kappa)/(2+kappa))

h = expand(132*(Phi12 - s*w + r))   # lcm(12,11,66,44) = 132
print(f"resultant(h, h', w) [deg 12 x 11] ... [{time.time()-t0:.0f}s]", flush=True)
D = resultant(h, diff(h, w), w)
print(f"  resultant done [{time.time()-t0:.0f}s]", flush=True)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D11 = sp.expand(D / cont)
P = sp.Poly(D11, s, r)
print("primitive D11: degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r), flush=True)
print("LA1 terms=76:", len(P.terms()) == 76, " degree=12:", P.total_degree() == 12, flush=True)
out["terms"] = len(P.terms()); out["degree"] = P.total_degree()
with open("atlas11_wall.txt","w") as f: f.write(str(D11))
supp = set(m for (m, _c) in P.terms())   # NB: Poly.terms() returns (monom, coeff) pairs
NN = 12
cone = [(i,j) for i in range(NN+1) for j in range(NN+1) if (NN-1)*i + NN*j <= NN*(NN-1)]
holes = [(i,j) for (i,j) in cone if (i,j) not in supp]
extra = sorted(supp - set(cone))
print("LA1 cone size:", len(cone), "[predict 79] | genuine holes:", holes, " [predict [(0,0),(1,0),(0,10)]]")
print("LA1 outside-cone support:", extra, " [predict NONE: (12,0) sits ON the geometric cone] | fictitious (0,12):", (0,12) not in supp)
out["holes"] = holes; out["extra"] = extra

z00 = int(D11.subs({s:0,r:0})); out["D11(0,0)"] = z00
print("LA2 D11(0,0) =", z00, "(predict 0)")
scr = sp.rem(D11.subs(r,0), s**2, s)
print("LA2 s^2 | D11(s,0):", sp.expand(scr) == 0, flush=True)
resN = int(sp.Poly(D, s, r).coeff_monomial(s**12))
law = R(132)**23 * 11**11 * R(1,12)**11
ratio = R(abs(resN), 1)/law
print(f"LA3 s^12 raw law ratio (==1): {ratio}", flush=True)
out["s12_raw_ratio"] = str(ratio)
s12 = int(P.coeff_monomial(s**12))
print("LA4 s^12 primitive =", s12, " sign (+):", s12 > 0, " factor:", sp.factorint(abs(s12)), flush=True)
out["s12"] = s12
content12 = sp.fraction(R(abs(resN), abs(s12)))[0] if s12 else None
sup12 = set(sp.factorint(content12)) if content12 > 1 else set()
print("LA4 content(12) =", content12, sp.factorint(content12) if content12 > 1 else "",
      " support subset {2,3,11}:", sup12 <= {2,3,11}, flush=True)
out["content12"] = str(content12)

fac = sp.factor(D11)
out["irreducible"] = bool(fac == D11 or fac == -D11)
print("LA1 irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)
tt = symbols("tt")
spara = p11.subs(w, tt); rpara = expand(tt*spara - Phi12.subs(w, tt))
out["param_on_wall"] = (sp.expand(D11.subs({s: spara, r: rpara})) == 0)
print("LA5 param identity:", out["param_on_wall"], flush=True)

p11p = diff(p11, w); p11pp = diff(p11p, w)
g12 = sp.gcd(p11p, p11pp)
print("LA6 gcd(p11',p11'') =", g12, " (1 => 10 simple ordinary cusps)", flush=True)
out["gcd_pp_ppp"] = str(g12)
# exact Sturm dance datum
stu = sp.sturm(sp.Poly(p11p, w), w)
nreal_sturm = sp.polys.polytools.count_roots(p11p, -sp.oo, sp.oo)
print(f"LA6 Sturm real roots of p11' over R: {nreal_sturm} [predict 2 - dance odd d]", flush=True)
out["sturm_real_cusps"] = int(nreal_sturm)

import mpmath as mp
mp.mp.dps = 110
def to_mpc(v):  return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
def synth_div(co, a):
    outq = [co[0]]
    for c in co[1:-1]:
        outq.append(c + a*outq[-1])
    return outq, co[-1] + a*outq[-1]

print("LA6 cusp contacts (predict 10, exactly 2 real, BOTH HIT - residual nonic 9/9 real):", flush=True)
cusps = []
for wv in sp.nroots(p11p, n=110, maxsteps=20000):
    sv = p11.subs(w, wv); rv = wv*sv - Phi12.subs(w, wv)
    cusps.append((wv, sv, rv))
out["real_cusps"] = 0; hit_report = []
for wv, sv, rv in cusps:
    real = abs(sp.im(wv)) < sp.Float("1e-90")
    if real: out["real_cusps"] += 1
    co = [to_mpc(c) for c in [-R(1,12), R(1,11), 0,0,0,0,0,0,0, -R(65,66), R(43,44), -sv, rv]]
    t_ = to_mpc(wv)
    g1, r1 = synth_div(co, t_); g2, r2 = synth_div(g1, t_); g, r3 = synth_div(g2, t_)
    resid = max(abs(r1), abs(r2), abs(r3))
    roots = mp.polyroots(g, maxsteps=6000, error=False)
    nreal = sum(1 for zz in roots if abs(mp.im(zz)) < mp.mpf("1e-25"))
    if real:
        hit_report.append((str(sp.N(wv,15)), str(sp.N(sv,15)), str(sp.N(rv,15)), nreal, mp.nstr(resid,3)))
        print(f"  REAL cusp t={sp.N(wv,12)} -> ({sp.N(sv,10)}, {sp.N(rv,10)}): residual nonic real roots {nreal}/9"
              f"  (audit {mp.nstr(resid,3)})  {'HIT' if nreal>=1 else 'NOT HIT'}", flush=True)
out["cusp_hits"] = hit_report
out["both_hit"] = all(x[3] >= 1 for x in hit_report) and len(hit_report) == 2
print("LA6 both real cusps HIT:", out["both_hit"],
      " | (sharpened 9/9 guess FALSIFIED - actual:", [x[3] for x in hit_report], "of 9)", flush=True)
coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(complex(sp.N(cusps[i][1]-cusps[j][1], 40))) < 1e-15 and abs(complex(sp.N(cusps[i][2]-cusps[j][2], 40))) < 1e-15)
print("cusp collisions:", coll, f"  [stage A {time.time()-t0:.0f}s]", flush=True)
out["cusp_collisions"] = int(coll); out["n_cusps"] = len(cusps)
json.dump(out, open("atlas11_stageA.json","w"), indent=1, default=str)
print("saved atlas11_stageA.json  GREEN:", len(P.terms())==76 and P.total_degree()==12 and z00==0
      and sp.expand(scr)==0 and ratio==1 and s12>0 and sup12<={2,3,11} and out["irreducible"]
      and out["param_on_wall"] and nreal_sturm==2 and out["both_hit"] and coll==0
      and set(holes)=={(0,0),(1,0),(0,10)} and extra==[] and len(cone)==79, flush=True)
