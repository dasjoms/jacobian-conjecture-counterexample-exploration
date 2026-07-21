"""
Note 9, stage A: wall of F7 (fiber OCTIC, d=7, even chamber - cone should return).
p7(w) = -w^7 + w^6 - (81/28)w^2 + (53/28)w
Predictions: wall deg 8; 6 cusps (roots of p7', expect 2 REAL); 15 nodes (1 crunode?);
budget 21; budget-balanced; S8 next stage.
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, factor, resultant, Rational as R, gcd as sgcd
import numpy as np

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p7 = -w**7 + w**6 - R(81,28)*w**2 + R(53,28)*w
Phi7 = expand(integrate(p7, w))
assert expand(p7 - (2*w - 3*w**2 + w*(1-w)*(w**5 - R(3,28)))) == 0
assert p7.subs(w,1) == -1 and Phi7.subs(w,1) == 0
kappa = sp.diff(p7,w).subs(w,1)
print(f"p7 = {p7}\nPhi7 = {Phi7}\nkappa = {kappa}, recipe a = {-(1+kappa)/(2+kappa)}", flush=True)
out["kappa"] = str(kappa)

h = expand(56*(Phi7 - s*w + r))   # lcm 8,7,28,56 = 56 -> integer scaling
print("resultant(h, h', w) [deg 8 x 7] ...", flush=True)
D = resultant(h, diff(h, w), w)
print(f"  raw resultant [{time.time()-t0:.0f}s]", flush=True)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D7 = sp.expand(D / cont)
P = sp.Poly(D7, s, r)
print("primitive D7: degree", P.total_degree(), "| terms:", len(P.terms()),
      "| deg s:", P.degree(s), "| deg r:", P.degree(r))
print("D7(0,0) =", D7.subs({s:0,r:0}))
out["D7"] = str(D7)
with open("atlas7_wall.txt","w") as f: f.write(str(D7))
print("D7 =", D7, flush=True)
fac = sp.factor(D7)
out["irreducible"] = bool(fac == D7 or fac == -D7)
print("irreducible:", out["irreducible"], f"[{time.time()-t0:.0f}s]", flush=True)

tt = symbols("tt")
spara = p7.subs(w, tt); rpara = expand(tt*spara - Phi7.subs(w, tt))
out["param_on_wall"] = (sp.expand(D7.subs({s: spara, r: rpara})) == 0)
print("param identity:", out["param_on_wall"], flush=True)
print("s^8 coeff factor:", sp.factorint(int(abs(P.coeff_monomial(s**8)))), flush=True)

p7p = diff(p7, w); p7pp = diff(p7p, w)
out["gcd(p7',p7'')"] = str(sp.Poly(sgcd(p7p, p7pp), w).as_expr())
print("gcd(p7', p7'') =", out["gcd(p7',p7'')"], flush=True)

cusps = []
for wv in sp.nroots(p7p, n=50):
    sv = p7.subs(w, wv); rv = wv*sv - Phi7.subs(w, wv)
    cusps.append((complex(sp.N(wv,50)), complex(sp.N(sv,50)), complex(sp.N(rv,50))))
print(f"cusp contacts: {len(cusps)}  [predict 6]")
for wv, sv, rv in cusps:
    real = abs(wv.imag) < 1e-30
    hc = [-R(1,8), R(1,7), 0, 0, 0, -R(27,28), R(53,56), -sv, rv]
    roots = np.roots([complex(c) for c in hc])
    others = [z for z in roots if abs(z - wv) > 1e-6]
    nreal = sum(1 for z in others if abs(z.imag) < 1e-9)
    print(f"  t={wv:.7f} -> ({sv:.7f},{rv:.7f})  {'REAL' if real else 'cplx'}  "
          f"residual real roots: {nreal}/5 {'=> missed over R' if real and nreal==0 else ('(hit, odd residuals)' if real else '')}", flush=True)
out["real_cusps"] = sum(1 for wv,_,_ in cusps if abs(wv.imag) < 1e-30)

coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
           if abs(cusps[i][1]-cusps[j][1]) < 1e-20 and abs(cusps[i][2]-cusps[j][2]) < 1e-20)
out["cusp_collisions"] = coll
print("cusp collisions:", coll, flush=True)
out["cusps"] = [[repr(a), repr(b), repr(c)] for a,b,c in cusps]
json.dump(out, open("atlas7_stageA.json","w"), indent=1)
print(f"stage A done [{time.time()-t0:.0f}s]", flush=True)
