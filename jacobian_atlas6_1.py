"""
Note 8, stage A: wall of F6 (fiber SEPTIC).
p6(w) = -w^6 + w^5 - (20/7)w^2 + (13/7)w   [d=6 explainer seed]
Phi6 = -w^7/7 + w^6/6 - (20/21)w^3 + (13/14)w^2
Predictions from note 7's pattern conjecture:
  wall degree 7, term count 27 = C(8,2)-1; 5 cusps (roots of p6'); 10 nodes;
  delta budget 15 = (7-1)(7-2)/2; no deeper strata (geometric emptiness check).
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, factor, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p6 = -w**6 + w**5 - R(20,7)*w**2 + R(13,7)*w
Phi6 = expand(integrate(p6, w))
assert expand(p6 - (2*w - 3*w**2 + w*(1-w)*(w**4 - R(1,7)))) == 0
assert p6.subs(w,1) == -1 and Phi6.subs(w,1) == 0
kappa = sp.diff(p6,w).subs(w,1)
print(f"p6 = {p6}\nPhi6 = {Phi6}\nkappa = {kappa}, recipe a = {-(1+kappa)/(2+kappa)}", flush=True)
out["kappa"] = str(kappa)

h = expand(210*(Phi6 - s*w + r))     # lcm(7,6,21,14) = 42 -> 210 safe integer scaling
print("computing resultant(h, h', w) [deg 7 x 6] ...", flush=True)
D = resultant(h, diff(h, w), w)
print(f"  raw resultant done [{time.time()-t0:.0f}s]", flush=True)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D6 = sp.expand(D / cont)
P = sp.Poly(D6, s, r)
print("primitive D6: total degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r),
      "| predicted 27 terms:", len(P.terms()) == 27)
print("D6(0,0) =", D6.subs({s:0,r:0}))
out["D6"] = str(D6)
with open("atlas6_wall.txt","w") as f: f.write(str(D6))
print("D6 =", D6, flush=True)
fac = sp.factor(D6)
out["irreducible"] = bool(fac == D6 or fac == -D6)
print("irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)

tt = symbols("tt")
spara = p6.subs(w, tt); rpara = expand(tt*spara - Phi6.subs(w, tt))
out["param_on_wall"] = (sp.expand(D6.subs({s: spara, r: rpara})) == 0)
print("parametrization on wall identically:", out["param_on_wall"], flush=True)

p6p = diff(p6, w); p6pp = diff(p6p, w)
out["gcd(p6',p6'')"] = str(sp.Poly(sgcd(p6p, p6pp), w).as_expr())
print("gcd(p6', p6'') =", out["gcd(p6',p6'')"], flush=True)

cusps = []
for wv in sp.nroots(p6p, n=50):
    sv = p6.subs(w, wv); rv = wv*sv - Phi6.subs(w, wv)
    cusps.append((complex(sp.N(wv,50)), complex(sp.N(sv,50)), complex(sp.N(rv,50))))
print(f"cusp contacts (roots of p6'): {len(cusps)}  [predict 5]", flush=True)
for wv, sv, rv in cusps:
    real = abs(wv.imag) < 1e-30
    # residual roots at the cusp: which real partners? (R-side whisker analysis)
    import numpy as np
    hc = [-R(1,7), R(1,6), 0, 0, -R(20,21), R(13,14), -sv, rv]
    roots = np.roots([float(c.m) if hasattr(c,'m') else complex(c) for c in hc])
    others = [z for z in roots if abs(z - wv) > 1e-6]
    nreal = sum(1 for z in others if abs(z.imag) < 1e-9)
    print(f"  t={wv:.7f} -> (s,r)=({sv:.7f},{rv:.7f})  {'REAL' if real else 'cplx'}  "
          f"residual real roots: {nreal}/4 {'=> missed over R!' if real and nreal==0 else ''}", flush=True)

coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(cusps[i][1]-cusps[j][1]) < 1e-20 and abs(cusps[i][2]-cusps[j][2]) < 1e-20)
out["cusp_collisions_(3,3)"] = coll
print("cusp collisions ((3,3) strata):", coll, flush=True)
out["cusps"] = [ [repr(a), repr(b), repr(c)] for a,b,c in cusps ]
json.dump(out, open("atlas6_stageA.json","w"), indent=1)
print(f"stage A done [{time.time()-t0:.0f}s]", flush=True)
