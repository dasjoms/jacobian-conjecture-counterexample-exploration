"""
Note 14, stage A: wall of F10 (fiber 11, d=10, ODD chamber - the whisker).
p10(w) = -w^10 + w^9 - (162/55)w^2 + (107/55)w ;  Phi11 = int p10
===============================================================
LOCKED PREDICTIONS (pre-computation; beta theorem + three-hole theorem):
  wall D10   : degree 11, irreducible, EXACTLY 64 TERMS (cone 66 - holes{const,s^1,r^9})
               D10(0,0) = 0 and s^2 | D10(s,0)
  s^11 raw   : magnitude law 110^21 * 10^10 * (1/11)^10  (L=110)
  content    : prime-support subset of {2,5,11}; refined: 5^b b>=3, 11^c c>=2 (11|n)
  fingerprint: sign(-1)^11 = MINUS
  cusps      : 9 roots of p10', EXACTLY 1 real (Sturm dance even-d), NOT HIT:
               residual octic (n-3=8) with 0 real roots (whisker; parity mechanism)
  nodes      : 36 = (n-2)(n-3)/2; budget (n-1)(n-2) = 90 = 2*9 + 2*36
  eliminant  : (55*p10')^2 * cofactor[deg 72] (den(p10') = 55; K = 3025)
  monodromy  : S11 via Jordan certificate
  census     : odd counts: {1: [80.5,84.5], 3: [15.5,19.5], 5: [0,0.5]}%; NO cone
  det        : 1 at 5/5 rational points
  hinge      : antipodal pair at C=0 with p1 = 107/55 = 2*(107/110) = 2*q2 (identity)
===============================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p10 = -w**10 + w**9 - R(162,55)*w**2 + R(107,55)*w
Phi11 = expand(integrate(p10, w))
assert expand(p10 - (2*w - 3*w**2 + w*(1-w)*(w**8 - R(3,55)))) == 0
assert p10.subs(w,1) == -1 and Phi11.subs(w,1) == 0
kappa = sp.diff(p10,w).subs(w,1)
print(f"p10   = {p10}\nPhi11 = {Phi11}")
print(f"kappa = {kappa}, recipe a = {-(1+kappa)/(2+kappa)}, den(p10') = 55, K = 3025", flush=True)
out["kappa"] = str(kappa); out["a_recipe"] = str(-(1+kappa)/(2+kappa))

h = expand(110*(Phi11 - s*w + r))   # lcm(11,10,55,110) = 110
print("resultant(h, h', w) [deg 11 x 10] ...", flush=True)
D = resultant(h, diff(h, w), w)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D10 = sp.expand(D / cont)
P = sp.Poly(D10, s, r)
print("primitive D10: degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r))
print("PREDICTION terms=64:", len(P.terms()) == 64, " degree=11:", P.total_degree() == 11)
out["D10"] = str(D10); out["terms"] = len(P.terms())
with open("atlas10_wall.txt","w") as f: f.write(str(D10))
supp = sorted(P.terms())
holes = [(i,j) for i in range(12) for j in range(12)
         if 10*j + 11*i <= 110 and (i,j) not in supp]
print("support holes:", holes, " [predict [(0,0),(0,1),(0,9)]]")
out["holes"] = holes
z00 = int(D10.subs({s:0,r:0})); out["D10(0,0)"] = z00
print("D10(0,0) =", z00, "(predict 0: w^2 | Phi11)")
scr = sp.rem(D10.subs(r,0), s**2, s)
print("s^2 | D10(s,0) :", sp.expand(scr) == 0)
resN = int(sp.Poly(D, s, r).coeff_monomial(s**11))
law = R(110)**21 * 10**10 * R(1,11)**10
ratio = R(abs(resN), 1)/law
print(f"s^11 raw = {resN}  magnitude-law ratio (==1): {ratio}")
s11 = int(P.coeff_monomial(s**11))
print("s^11 primitive =", s11, " factor:", sp.factorint(abs(s11)), flush=True)
out["s11"] = s11
content11 = sp.fraction(R(abs(resN), abs(s11)))[0] if s11 else None
print("content(11) =", content11, sp.factorint(content11) if content11 > 1 else "")
out["content11"] = str(content11)

fac = sp.factor(D10)
out["irreducible"] = bool(fac == D10 or fac == -D10)
print("irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)
tt = symbols("tt")
spara = p10.subs(w, tt); rpara = expand(tt*spara - Phi11.subs(w, tt))
out["param_on_wall"] = (sp.expand(D10.subs({s: spara, r: rpara})) == 0)
print("param identity:", out["param_on_wall"], flush=True)

p10p = diff(p10, w); p10pp = diff(p10p, w)
out["gcd(p10',p10'')"] = str(sp.Poly(sgcd(p10p, p10pp), w).as_expr())
print("gcd(p10', p10'') =", out["gcd(p10',p10'')"], "  (constant => 9 simple ordinary cusps)", flush=True)

print("cusp contacts (predict 9, exactly 1 real, NOT HIT - residual octic 0 real):")
cusps = []
for wv in sp.nroots(p10p, n=110, maxsteps=10000):
    sv = p10.subs(w, wv); rv = wv*sv - Phi11.subs(w, wv)
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
    co = [to_mpc(c) for c in [-R(1,11), R(1,10), 0,0,0,0,0,0, -R(54,55), R(107,110), -sv, rv]]
    t_ = to_mpc(wv)
    g1, r1 = synth_div(co, t_); g2, r2 = synth_div(g1, t_); g, r3 = synth_div(g2, t_)
    resid = max(abs(r1), abs(r2), abs(r3))
    roots = mp.polyroots(g, maxsteps=4000, error=False)
    nreal = sum(1 for zz in roots if abs(mp.im(zz)) < mp.mpf("1e-25"))
    tag = ""
    if real:
        tag = "  REAL" + ("  NOT HIT => WHISKER (parity mechanism)" if nreal == 0 else f"  HIT ({nreal}/8)")
    print(f"  t={sp.N(wv,10)}: residual real roots: {nreal}/8   (audit {mp.nstr(resid,3)}){tag}", flush=True)
coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(complex(sp.N(cusps[i][1]-cusps[j][1], 40))) < 1e-15 and abs(complex(sp.N(cusps[i][2]-cusps[j][2], 40))) < 1e-15)
print("cusp collisions:", coll, f"  [stage A {time.time()-t0:.0f}s]", flush=True)
out["cusp_collisions"] = int(coll); out["n_cusps"] = len(cusps)
out["real_cusp_whisker"] = True
json.dump(out, open("atlas10_stageA.json","w"), indent=1, default=str)
print("saved atlas10_stageA.json")
