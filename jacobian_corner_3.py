"""
NOTE 15, stage 2d: EXACT count-window renders + node census reconciliation.
No new locks: pure cartography (audit-grade) + census parse + d=4 r-solve.
Windows:
 d=2: [-2.5,0.6] x [-2.5,0.6]   (cusp-interior missed region, tip (1/3,1/27))
 d=3: [-1.6,-0.4]^2             (crunode EXACT (-1,-1))
 d=4: [-1.7,-0.2] x [-1.7,0.1]  (whisker-hunt; myth (-4/3,-4/3) burial plot)
 d=4b:[-0.1,+0.5] x [-0.2,+0.3] (ghost cusp region)
 d=5: [-1.15,-0.35]^2           (cone vertex region)
"""
import json, time
import numpy as np
import sympy as sp
from sympy import symbols, expand, diff, integrate, resultant, Rational as R, Poly, factor

w, s, r = symbols("w s r")
t0 = time.time()
out = {}

def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

def exact_count(d, sv, rv):
    p = seed(d); Phi = integrate(p, w)
    hd = Poly(sp.expand((Phi - w*sv + rv)), w)
    return len(hd.real_roots())

def render(d, x0, x1, y0, y1, n, tag):
    S = [x0 + (x1-x0)*i/(n-1) for i in range(n)]
    Rv = [y0 + (y1-y0)*i/(n-1) for i in range(n)]
    G = np.zeros((n, n), dtype=int)
    for i, rv in enumerate(Rv):
        for j, sv in enumerate(S):
            G[i, j] = exact_count(d, R(f"{sv:.9f}"), R(f"{rv:.9f}"))
    np.savez(f"jcorner_map_{tag}.npz", S=np.array(S), R=np.array(Rv), G=G)
    uniq, cnt = np.unique(G, return_counts=True)
    print(f"  [{tag}] d={d} counts {dict(zip(uniq.tolist(), cnt.tolist()))}", flush=True)
    return G, S, Rv

print("="*84); print("exact count windows"); print("-"*84)
G2,  S2,  R2  = render(2, -2.5,  0.6, -2.5,  0.6, 40, "d2")
G3,  S3,  R3  = render(3, -1.6, -0.4, -1.6, -0.4, 44, "d3")
G4,  S4,  R4  = render(4, -1.7, -0.2, -1.7,  0.1, 44, "d4")
G4b, S4b, R4b = render(4, -0.10, 0.50, -0.20, 0.30, 44, "d4b")
G5,  S5,  R5  = render(5, -1.15, -0.35, -1.15, 0.15, 44, "d5")

# ---- d=5 vertex from the map: leftmost (in s) count-0 cell ----
i0 = np.argwhere(G5 == 0)
if len(i0):
    js = [j for i, j in i0]
    jmin = min(js)
    print(f"  d=5: count-0 region leftmost column s ~ {S5[jmin]:.4f}")
# diagonal trace d=5: count along r=s (nearest cells)
tr = [(sv, G5[i, np.argmin([abs(x-sv) for x in S5])]) for i, sv in enumerate(R5) if abs(sv - S5[np.argmin([abs(x-sv) for x in S5])]) < 0.01]
print(f"  d=5 diagonal counts (r=s cells): {sorted(set(c for _, c in tr))}")

# ---- d=4 r-solve at the real node-s ----
print("="*84); print("d=4 ghost-node r-solve + whisker diagnosis"); print("-"*84)
p4 = seed(4); Phi4 = integrate(p4, w)
h4 = expand((Phi4 - s*w + r)*10)
D4 = resultant(h4, diff(h4, w), w)
D4 = expand(D4/sp.gcd_list(Poly(D4, s, r).coeffs()))
node_cub = Poly(40960000*w**3 - 112883200*w**2 + 105406221*w - 33445465, w)
# same cubic in s
nc = Poly(40960000*s**3 - 112883200*s**2 + 105406221*s - 33445465, s)
sreal = sp.N(nc.real_roots()[0], 40)
rvals = sp.Poly(D4.subs(s, sp.Float(sreal, 40)), r).nroots(30)
grad = [rv for rv in rvals if abs(complex(sp.N(diff(D4, s).subs({s: sreal, r: rv}), 25))) < 1e-10]
print(f"  node-s = {sreal}; r candidates with dD/ds ~ 0: {[str(z)[:24] for z in grad]}")
out["d4_ghost_node"] = {"s": str(sreal), "r": str(grad[0]) if grad else None}

# whisker-root candidates for d=4 -- from the count maps:
# whisker = arc bounding count-drop pair; find for d=4 the 1-real region's extent
def region_stats(G, S, Rv, val, tag):
    ij = np.argwhere(G == val)
    if not len(ij):
        print(f"  [{tag}] count-{val} region: EMPTY"); return None
    rs = sorted(set(i for i, j in ij)); cs = sorted(set(j for i, j in ij))
    st = {"s_min": S[min(cs)], "s_max": S[max(cs)], "r_min": Rv[min(rs)], "r_max": Rv[max(rs)]}
    print(f"  [{tag}] count-{val} bbox s:[{st['s_min']:.3f},{st['s_max']:.3f}] r:[{st['r_min']:.3f},{st['r_max']:.3f}]", flush=True)
    return st
region_stats(G4, S4, R4, 1, "d4-corner-win")     # 1-real (whisker-adjacent) region
region_stats(G4b, S4b, R4b, 1, "d4-ghost-win")
region_stats(G4, S4, R4, 5, "d4-busy")

json.dump(out, open("jcorner_stage2d.json", "w"), indent=1)
print(f"saved jcorner_stage2d.json + maps [{time.time()-t0:.0f}s]")
