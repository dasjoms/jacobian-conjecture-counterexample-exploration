"""
Note 7, stage A: wall of F5 (fiber SEXTIC).
p5(w) = 2w - 3w^2 + w(1-w)(w^3 - 1/5)  (d=5 explainer seed)
Phi5 = int p5 ; h(w;s,r) = Phi5 - s w + r ; wall D5 = resultant(h, h', w).
Genus bookkeeping target: rational sextic => sum of delta-invariants == (6-1)(6-2)/2 == 10.
Note 4 emptiness theorem predicts: no (2,2,2) [triple point], no (3,3), no (4,2), no (6).
"""
import sympy as sp, json, time
from sympy import symbols, diff, integrate, expand, factor, resultant, Rational as R, gcd as sgcd

t0 = time.time()
w, s, r = symbols("w s r")
out = {}

p5 = -w**5 + w**4 - R(14,5)*w**2 + R(9,5)*w
Phi5 = expand(integrate(p5, w))
assert expand(p5 - (2*w - 3*w**2 + w*(1-w)*(w**3 - R(1,5)))) == 0
assert p5.subs(w,1) == -1 and Phi5.subs(w,1) == 0
kappa = sp.diff(p5,w).subs(w,1)
out["kappa"] = str(kappa); out["a_recipe"] = str(-(1+kappa)/(2+kappa))
print(f"p5 = {p5}\nPhi5 = {Phi5}\nkappa = {kappa}, a = {-(1+kappa)/(2+kappa)}", flush=True)

# ---- wall discriminant: clear denominators first (30 = lcm) ----
h = expand(30*(Phi5 - s*w + r))          # integer coeffs now
print("computing resultant(h, h', w) [deg 6 x 5] ...", flush=True)
D = resultant(h, diff(h, w), w)
print(f"  raw resultant done [{time.time()-t0:.0f}s]", flush=True)
cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
D5 = sp.expand(D / cont)
print("primitive D5 degree:", sp.Poly(D5, s, r).total_degree(),
      " terms:", len(sp.Poly(D5, s, r).terms()),
      " deg in s:", sp.Poly(D5, s, r).degree(s), " deg in r:", sp.Poly(D5, s, r).degree(r))
print("D5(0,0) =", D5.subs({s:0, r:0}), "-> C=0 in wall")
out["D5"] = str(D5)
with open("atlas5_wall.txt","w") as f: f.write(str(D5))
print("D5 =", D5, flush=True)
print(f"[{time.time()-t0:.0f}s]", flush=True)
try:
    fac = sp.factor(D5)
    out["irreducible"] = bool(fac == D5 or fac == -D5)
    print("irreducible over Q(s,r):", out["irreducible"], flush=True)
except Exception as e:
    out["irreducible"] = f"factor raised {e}"
print(f"[{time.time()-t0:.0f}s]", flush=True)

# ---- parametrization check: t |-> (s(t), r(t)) satisfies D5 == 0 identically ----
tt = symbols("tt")
spara = p5.subs(w, tt)
rpara = expand(tt*spara - Phi5.subs(w, tt))
chk = sp.expand(D5.subs({s: spara, r: rpara}))
out["param_on_wall"] = (chk == 0)
print("parametrization satisfies D5 identically:", chk == 0, flush=True)

# ---- cusp contacts: roots of p5' ----
p5p = diff(p5, w); p5pp = diff(p5p, w)
gcusp = sp.Poly(sgcd(p5p, p5pp), w)
out["gcd(p5',p5'')"] = str(gcusp.as_expr())
print("gcd(p5', p5'') =", gcusp.as_expr(), " -> >=4-fold contacts:", gcusp.degree() > 0, flush=True)

cusps = []
for wv in sp.nroots(p5p, n=50):
    sv = p5.subs(w, wv); rv = wv*sv - Phi5.subs(w, wv)
    cusps.append((complex(sp.N(wv,50)), complex(sp.N(sv,50)), complex(sp.N(rv,50))))
print(f"cusp contacts (roots of p5'): {len(cusps)}")
for wv, sv, rv in cusps:
    real = abs(wv.imag) < 1e-30
    print(f"  t={wv} -> (s,r)=({sv},{rv})  {'REAL' if real else 'complex'}", flush=True)

# cusp collisions = (3,3) strata
coll = []
for i in range(len(cusps)):
    for j in range(i+1, len(cusps)):
        if abs(cusps[i][1]-cusps[j][1]) < 1e-20 and abs(cusps[i][2]-cusps[j][2]) < 1e-20:
            coll.append((i,j))
out["cusp_collisions_(3,3)"] = len(coll)
print("cusp-point collisions ((3,3) strata):", len(coll), flush=True)

# ---- (3,2,1): at each cusp t0, residual cubic h/(w-t0)^3 and its discriminant ----
print("--- (3,2,1) tacnodes: cusp + extra double root ---")
tac = []
h_sym = sp.expand(30*(Phi5 - s*w + r))
for wv, sv, rv in cusps:
    hc = sp.Poly(h_sym.subs({s: sv, r: rv}), w)
    # quotient by (w - wv)^3 numerically via mpmath-friendly: use exact conjugates?
    # numeric route: divide polynomial coefficients
    import numpy as np
    coeffs = [complex(sp.N(c,50)) for c in hc.all_coeffs()]
    roots = np.roots(coeffs)
    others = [z for z in roots if abs(z - wv) > 1e-8]
    # find duplicate among 'others'
    dup = None
    for i in range(len(others)):
        for j in range(i+1, len(others)):
            if abs(others[i]-others[j]) < 1e-7:
                dup = others[i]
    mult0 = sum(1 for z in roots if abs(z - wv) <= 1e-6)
    print(f"  cusp t={wv:.6f}: root mult at t0 = {mult0}, others = {[complex(round(z.real,5), round(z.imag,5)) for z in others]}"
          + (f"  EXTRA DOUBLE at {dup:.6f} -> (3,2,1)!" if dup is not None else ""), flush=True)
    if dup is not None:
        tac.append((wv, sv, rv, dup))
out["tacnodes_(3,2,1)"] = len(tac)

with open("atlas5_stageA.json","w") as f:
    json.dump(out, f, indent=1)
print(f"stage A core done [{time.time()-t0:.0f}s]", flush=True)
