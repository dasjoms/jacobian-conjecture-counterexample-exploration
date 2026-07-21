"""
Note 11, stage A: wall of F9 (fiber DECIC, d=9, EVEN chamber - the cone returns).
p9(w) = -w^9 + w^8 - (44/15)w^2 + (29/15)w ;  Phi10 = int p9
===============================================================
LOCKED PREDICTIONS (pre-computation; beta(10)=0 certified by holes_0):
  wall D9    : degree 10, irreducible, EXACTLY 53 TERMS (cone 56 - holes{const,s^1,r^8})
               D9(0,0) = 0 and s^2 | D9(s,0) (holes 1+2: w^2 | Phi mechanism)
  s^10 raw   : magnitude law 2^10 * 3^56 * 5^10 (L=90); content unknown (open law)
  cusps      : 8 roots of p9', EXACTLY 2 real (Sturm dance), BOTH HIT
               (residual septic: odd degree => >= 1 real root ALWAYS - mechanism)
  nodes      : 28 = (n-2)(n-3)/2; budget 8+28 = 36; 1 crunode (dance), 3 acnodes (staircase)
  eliminant  : (15*p9')^2 * cofactor[deg 56] (den(p9') = 15)
  monodromy  : S10 via Jordan certificate (transposition + transitivity on pairs)
  census     : even counts, missed CONE ~ 8.7% (echo of 8.69/8.74)
  det        : 1 at 5/5 rational points (and retro d=6, d=7)
===============================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p9 = -w**9 + w**8 - R(44,15)*w**2 + R(29,15)*w
Phi10 = expand(integrate(p9, w))
assert expand(p9 - (2*w - 3*w**2 + w*(1-w)*(w**7 - R(1,15)))) == 0
assert p9.subs(w,1) == -1 and Phi10.subs(w,1) == 0
kappa = sp.diff(p9,w).subs(w,1)
print(f"p9    = {p9}\nPhi10 = {Phi10}")
print(f"kappa = {kappa}, recipe a = {-(1+kappa)/(2+kappa)}, den(p9') = 15", flush=True)
out["kappa"] = str(kappa)

h = expand(90*(Phi10 - s*w + r))   # lcm(10,9,45,30) = 90
print("resultant(h, h', w) [deg 10 x 9] ...", flush=True)
D = resultant(h, diff(h, w), w)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D9 = sp.expand(D / cont)
P = sp.Poly(D9, s, r)
print("primitive D9: degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r))
print("PREDICTION terms=53:", len(P.terms()) == 53, " degree=10:", P.total_degree() == 10)
out["D9"] = str(D9); out["terms"] = len(P.terms())
with open("atlas9_wall.txt","w") as f: f.write(str(D9))
z00 = int(D9.subs({s:0,r:0})); out["D9(0,0)"] = z00
print("D9(0,0) =", z00, "(predict 0: w^2 | Phi10)")
scr = sp.rem(D9, s**2, s)
print("s^2 | D9(s,0) :", sp.expand(scr.subs(r,0)) == 0)
resN = int(sp.Poly(D, s, r).coeff_monomial(s**10))
an10 = -R(1,10)
law = R(90)**19 * 9**9 * abs(an10)**9
ratio = R(abs(resN), 1)/law
print(f"s^10 raw = {resN}  magnitude-law ratio (==1): {ratio}")
s10 = int(P.coeff_monomial(s**10))
print("s^10 primitive =", s10, " factor:", sp.factorint(abs(s10)), flush=True)
out["s10"] = s10
content10 = sp.fraction(R(abs(resN), abs(s10)))[0] if s10 else None
print("content(10) =", content10, sp.factorint(content10) if content10 > 1 else "")
out["content10"] = str(content10)

fac = sp.factor(D9)
out["irreducible"] = bool(fac == D9 or fac == -D9)
print("irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)
tt = symbols("tt")
spara = p9.subs(w, tt); rpara = expand(tt*spara - Phi10.subs(w, tt))
out["param_on_wall"] = (sp.expand(D9.subs({s: spara, r: rpara})) == 0)
print("param identity:", out["param_on_wall"], flush=True)

p9p = diff(p9, w); p9pp = diff(p9p, w)
out["gcd(p9',p9'')"] = str(sp.Poly(sgcd(p9p, p9pp), w).as_expr())
print("gcd(p9', p9'') =", out["gcd(p9',p9'')"], "  (constant => 8 simple ordinary cusps)", flush=True)

print("cusp contacts (predict 8, exactly 2 real, BOTH HIT - residual septic is ODD degree):")
cusps = []
for wv in sp.nroots(p9p, n=110, maxsteps=8000):
    sv = p9.subs(w, wv); rv = wv*sv - Phi10.subs(w, wv)
    cusps.append((wv, sv, rv))

import mpmath as mp
mp.mp.dps = 110
def to_mpc(v):  return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
def synth_div(co, a):
    outq = [co[0]]
    for c in co[1:-1]:
        outq.append(c + a*outq[-1])
    return outq, co[-1] + a*outq[-1]
out["real_cusps"] = 0
for wv, sv, rv in cusps:
    real = abs(sp.im(wv)) < sp.Float("1e-90")
    if real: out["real_cusps"] += 1
    co = [to_mpc(c) for c in [-R(1,10), R(1,9), 0,0,0,0,0, -R(44,45), R(29,30), -sv, rv]]
    t_ = to_mpc(wv)
    g1, r1 = synth_div(co, t_); g2, r2 = synth_div(g1, t_); g, r3 = synth_div(g2, t_)
    resid = max(abs(r1), abs(r2), abs(r3))
    roots = mp.polyroots(g, maxsteps=3000, error=False)
    nreal = sum(1 for zz in roots if abs(mp.im(zz)) < mp.mpf("1e-25"))
    tag = ""
    if real:
        tag = "  REAL" + ("  => WHISKER?!" if nreal == 0 else f"  HIT ({nreal}/7 real residual, septic odd => >=1 automatic)")
    print(f"  t={sp.N(wv,10)}: residual real roots: {nreal}/7   (audit {mp.nstr(resid,3)}){tag}", flush=True)
coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(complex(sp.N(cusps[i][1]-cusps[j][1], 40))) < 1e-15 and abs(complex(sp.N(cusps[i][2]-cusps[j][2], 40))) < 1e-15)
print("cusp collisions:", coll, f"  [stage A {time.time()-t0:.0f}s]", flush=True)
out["cusp_collisions"] = int(coll)
out["n_cusps"] = len(cusps)
json.dump(out, open("atlas9_stageA.json","w"), indent=1, default=str)
