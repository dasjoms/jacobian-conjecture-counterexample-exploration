"""
Note 10, stage A: wall of F8 (fiber NONIC, d=8, ODD chamber - whisker protocol).
p8(w) = -w^8 + w^7 - (35/12)w^2 + (23/12)w ;  Phi9 = int p8
===============================================================
LOCKED PREDICTIONS (before any computation tonight):
  wall D8    : degree 9, EXACTLY 42 terms (6,10,14,20,26,34 -> +8), irreducible
  cusps      : 7 roots of p8', EXACTLY 1 real (dance 1,2,1,2,1,2 -> 1),
               residual sextic all-complex => WHISKER (odd chamber)
  nodes      : 21 = (n-2)(n-3)/2, budget 28 = 7 + 21 balanced
               0 crunodes (crunode dance 0,1,0,1 -> 0); acnodes ~ 3 (guess)
  eliminant  : (12*p8')^2 * cofactor[deg 42]  (denominator of p8' = 12)
  monodromy  : S9 (|G| = 362880)
  census     : odd counts only {1,3,5,7,9}; open-missed 0.0000% (whisker tip only)
  rates      : fold 1/2, cusp 2/3 (5th chamber each)
===============================================================
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p8 = -w**8 + w**7 - R(35,12)*w**2 + R(23,12)*w
Phi9 = expand(integrate(p8, w))
assert expand(p8 - (2*w - 3*w**2 + w*(1-w)*(w**6 - R(1,12)))) == 0
assert p8.subs(w,1) == -1 and Phi9.subs(w,1) == 0
kappa = sp.diff(p8,w).subs(w,1)
print(f"p8   = {p8}\nPhi9 = {Phi9}")
print(f"kappa = {kappa}, recipe a = {-(1+kappa)/(2+kappa)}, den(p8') = 12", flush=True)
out["kappa"] = str(kappa)

h = expand(72*(Phi9 - s*w + r))   # lcm(9,8,36,24) = 72 -> integer scaling
print("resultant(h, h', w) [deg 9 x 8] ...", flush=True)
D = resultant(h, diff(h, w), w)
print(f"  raw resultant [{time.time()-t0:.1f}s]", flush=True)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D8 = sp.expand(D / cont)
P = sp.Poly(D8, s, r)
print("primitive D8: degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r))
print("PREDICTION terms=42:", len(P.terms()) == 42, " degree=9:", P.total_degree() == 9)
out["D8"] = str(D8); out["terms"] = len(P.terms())
with open("atlas8_wall.txt","w") as f: f.write(str(D8))
print("D8 =", D8, flush=True)
z00 = int(D8.subs({s:0,r:0})); out["D8(0,0)"] = z00
print("D8(0,0) =", z00, " factor:", sp.factorint(abs(z00)), flush=True)
s9 = int(P.coeff_monomial(s**9)); out["s9"] = s9
print("s^9 coeff =", s9, " factor:", sp.factorint(abs(s9)), flush=True)

fac = sp.factor(D8)
out["irreducible"] = bool(fac == D8 or fac == -D8)
print("irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)

tt = symbols("tt")
spara = p8.subs(w, tt); rpara = expand(tt*spara - Phi9.subs(w, tt))
out["param_on_wall"] = (sp.expand(D8.subs({s: spara, r: rpara})) == 0)
print("param identity:", out["param_on_wall"], flush=True)

p8p = diff(p8, w); p8pp = diff(p8p, w)
out["gcd(p8',p8'')"] = str(sp.Poly(sgcd(p8p, p8pp), w).as_expr())
print("gcd(p8', p8'') =", out["gcd(p8',p8'')"], "  (1 => 7 simple ordinary cusps)", flush=True)

# --- cusp census at high precision (np.roots-free, audit by residuals) ---
print("cusp contacts (predict 7, exactly 1 real):")
cusps = []
for wv in sp.nroots(p8p, n=110, maxsteps=4000):
    sv = p8.subs(w, wv); rv = wv*sv - Phi9.subs(w, wv)
    cusps.append((wv, sv, rv))
    print(f"  t={sp.N(wv,12)}  ->  ({sp.N(sv,12)}, {sp.N(rv,12)})", flush=True)

import mpmath as mp
mp.mp.dps = 110
def to_mpc(v):  return mp.mpc(mp.mpf(str(sp.re(v))), mp.mpf(str(sp.im(v))))
def synth_div(co, a):   # divide poly(co, descending) by (w - a); return (quotient, remainder)
    outq = [co[0]]
    for c in co[1:-1]:
        outq.append(c + a*outq[-1])
    return outq, co[-1] + a*outq[-1]
out["real_cusps"] = 0
for wv, sv, rv in cusps:
    real = abs(sp.im(wv)) < sp.Float("1e-90")
    if real: out["real_cusps"] += 1
    co = [to_mpc(c) for c in [-R(1,9), R(1,8), 0,0,0,0, -R(35,36), R(23,24), -sv, rv]]
    t_ = to_mpc(wv)
    g1, r1 = synth_div(co, t_)
    g2, r2 = synth_div(g1, t_)
    g, r3 = synth_div(g2, t_)          # residual sextic h/(w-t)^3; ri should be ~1e-100
    resid = max(abs(r1), abs(r2), abs(r3))
    roots = mp.polyroots(g, maxsteps=2000, error=False)
    nreal = sum(1 for z in roots if abs(mp.im(z)) < mp.mpf("1e-25"))
    tag = ""
    if real: tag = "  REAL" + ("  => WHISKER (missed over R)" if nreal == 0 else f"  (hit, {nreal} real residual)")
    print(f"  t={sp.N(wv,10)}: residual real roots: {nreal}/6   (synth-div audit {mp.nstr(resid,3)}){tag}", flush=True)
    out[f"cusp_{sp.N(sp.re(wv),10)}"] = {"real": bool(real), "resid_real": nreal}

coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(complex(sp.N(cusps[i][1]-cusps[j][1], 40))) < 1e-15 and abs(complex(sp.N(cusps[i][2]-cusps[j][2], 40))) < 1e-15)
out["cusp_collisions"] = int(coll)
print("cusp collisions (|.| metric, distinct cusp points):", coll, f"  [stage A {time.time()-t0:.0f}s]", flush=True)
json.dump(out, open("atlas8_stageA.json","w"), indent=1, default=str)
