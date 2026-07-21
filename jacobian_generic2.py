"""
Note 10, stage 2:
(B) GENERIC-SEED ATLAS: 8 normalized deg-4 seeds (random (c3,c4) -> p1, c2 solved)
    through wall/strata/budget/emptiness; S5 monodromy + slopes on one fresh seed;
    d=3 chamber S4 monodromy.
(C) FORCED SWALLOWTAIL: deg-5 seed with p'(1/2) = p''(1/2) = 0 -> genuine A3
    (swallowtail) stratum: fiber root of multiplicity 4; escape slope ~ 3/4?
"""
import sympy as sp, json
from sympy import symbols, diff, integrate, expand, resultant, Rational as R, gcd as sgcd
import numpy as np

w, s, r = symbols("w s r")

def wall_of(p):
    Phi = expand(integrate(p, w))
    coeffs = sp.Poly(expand(Phi - s*w + r), w).all_coeffs()
    L = 1
    for c in coeffs:
        L = sp.ilcm(int(L), int(sp.denom(c)))
    h = expand(L*(Phi - s*w + r))
    D = resultant(h, diff(h, w), w)
    cont = sp.gcd_list([c for c in sp.Poly(D, s, r).coeffs()])
    return sp.expand(D/cont), Phi

def strata_of(p, Phi, prec=60):
    cusps = []
    for wv in sp.nroots(diff(p, w), n=50):
        sv = p.subs(w, wv); rv = wv*sv - Phi.subs(w, wv)
        cusps.append((complex(sp.N(wv,50)), complex(sp.N(sv,50)), complex(sp.N(rv,50))))
    w1, w2 = symbols("w1 w2")
    eq1 = sp.expand(sp.cancel((p.subs(w, w2) - p.subs(w, w1))/(w2 - w1)))
    eq2 = sp.expand(sp.cancel((Phi.subs(w, w2) - Phi.subs(w, w1))/(w2 - w1) - p.subs(w, w1)))
    gb = sp.groebner([eq1, eq2], w1, w2, order="lex")
    elim = None
    for g in gb.polys:
        if set(g.free_symbols) == {w2}:
            elim = g.expr
    cs, elimsf = sp.Poly(elim, w2).sqf_list()
    bit_elim = None
    for fac, mult in elimsf:
        if mult == 1 and (bit_elim is None or sp.degree(fac) > sp.degree(bit_elim)):
            bit_elim = fac
    nodes, gaps = [], []
    if bit_elim is not None and sp.degree(bit_elim) > 0:
        roots = sp.nroots(sp.Poly(bit_elim, w2), n=prec, maxsteps=3000)
        imgs = []
        for cv in roots:
            sv = complex(sp.N(p.subs(w, cv), prec)); rv = complex(sp.N(cv*p.subs(w, cv) - Phi.subs(w, cv), prec))
            imgs.append((complex(sp.N(cv, prec)), sv, rv))
        used = set()
        for i in range(len(imgs)):
            if i in used: continue
            best, bj = 1e30, -1
            for j in range(i+1, len(imgs)):
                if j in used: continue
                d = abs(imgs[i][1]-imgs[j][1]) + abs(imgs[i][2]-imgs[j][2])
                if d < best: best, bj = d, j
            used.add(i); used.add(bj)
            nodes.append((imgs[i], imgs[bj])); gaps.append(best)
    coll = sum(1 for i in range(len(cusps)) for j in range(i+1, len(cusps))
               if abs(cusps[i][1]-cusps[j][1]) < 1e-20 and abs(cusps[i][2]-cusps[j][2]) < 1e-20)
    trip = 0
    images = [(a[1], a[2]) for a, b in nodes]
    for i in range(len(images)):
        for j in range(i+1, len(images)):
            if abs(images[i][0]-images[j][0]) < 1e-15 and abs(images[i][1]-images[j][1]) < 1e-15:
                trip += 1
    overl = sum(1 for (a,b) in nodes for csv,c_rv in [(c[1],c[2]) for c in cusps]
                if abs(a[1]-csv) < 1e-15 and abs(a[2]-c_rv) < 1e-15)
    g = sgcd(diff(p, w), diff(p, w, 2))
    return cusps, nodes, gaps, coll, trip, overl, sp.degree(g, w)

print("=== (B) generic deg-4 seed atlas (fiber 5 each) ===")
draws = {"tower (1,-1)": (1,-1), "(2,-1)": (2,-1), "(-1,1)": (-1,1), "(2,2)": (2,2),
         "(-2,-1)": (-2,-1), "(1,2)": (1,2), "(3,-2)": (3,-2), "(-3,2)": (-3,2)}
report = {}
for name, (c3, c4) in draws.items():
    c3, c4 = R(c3), R(c4)
    p1 = 2 + c3/2 + 4*c4/5
    c2 = -3 - R(3,2)*c3 - R(9,5)*c4
    p = p1*w + c2*w**2 + c3*w**3 + c4*w**4
    norm_ok = (p.subs(w, 1) == -1) and (integrate(p, w).subs(w, 1) == 0)
    kap = diff(p, w).subs(w, 1)
    if not norm_ok or kap == -2:
        print(f"  {name}: DEGENERATE normalization/kappa"); continue
    D, Phi = wall_of(p)
    Ps = sp.Poly(D, s, r)
    cusps, nodes, gaps, coll, trip, overl, deepgcd = strata_of(p, Phi)
    fbdeg = 5
    budget_target = (fbdeg-1)*(fbdeg-2)//2
    n_real_cusps = sum(1 for cv,_,_ in cusps if abs(cv.imag) < 1e-25)
    real_nodes = [pr for pr in nodes if abs(pr[0][1].imag) < 1e-15 and abs(pr[0][2].imag) < 1e-15]
    balanced = (len(cusps) + len(nodes) + 2*coll + 2*trip + 2*overl == budget_target)
    empty = (coll == 0 and trip == 0 and overl == 0 and deepgcd == 0)
    report[name] = {"p": str(p), "wall_deg": Ps.total_degree(), "wall_terms": len(Ps.terms()),
                    "cusps": len(cusps), "nodes": len(nodes), "real_cusps": n_real_cusps,
                    "real_nodes": len(real_nodes), "balanced": balanced, "empty": empty,
                    "max_gap": float(max(gaps)) if gaps else None}
    msg = (f"  {name:16s}: wall deg {Ps.total_degree()} terms {len(Ps.terms())} | cusps {len(cusps)} "
           f"({n_real_cusps}R) nodes {len(nodes)} ({len(real_nodes)}R) | budget {len(cusps)}+{len(nodes)} "
           f"vs {budget_target} {'OK' if balanced else 'FAIL'} | empty {'OK' if empty else 'FAIL c%d/t%d/o%d/g%d' % (coll,trip,overl,deepgcd)}")
    if not nodes:
        msg += " | NO NODES?!"
    print(msg, flush=True)
json.dump(report, open("generic_atlas.json","w"), indent=1)

print("\n=== (B2) S5 spot-check on fresh generic seed (c3,c4)=(2,-1) ===")
c3, c4 = R(2), R(-1)
p1 = 2 + c3/2 + 4*c4/5; c2 = -3 - R(3,2)*c3 - R(9,5)*c4
pG = p1*w + c2*w**2 + c3*w**3 + c4*w**4
DG, PhiG = wall_of(pG)
DGl = sp.lambdify((s, r), DG, "numpy")
def hG(sv, rv):
    return [complex(c) for c in sp.Poly(expand((PhiG - s*w + r).subs({s: sv, r: rv})), w).all_coeffs()]
pGn = sp.lambdify(w, pG, "numpy")
sf = float(pG.subs(w, R(3,5))); rf = float(R(3,5)*pG.subs(w, R(3,5)) - PhiG.subs(w, R(3,5)))
loops = {"fold": lambda t: (sf + 0.005*np.exp(1j*t)/np.sqrt(2), rf + 1j*0.005*np.exp(1j*t)/np.sqrt(2)),
         "s=30e^{it},r=1": lambda t: (30*np.exp(1j*t), 1.0),
         "r=30e^{it},s=0.2": lambda t: (0.2, 30*np.exp(1j*t))}
N = 5
def track(path, steps):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    md = min(abs(DGl(a, b)) for a, b in pts)
    def greedy(prev, nxt):
        used = [False]*N; out = [None]*N
        for a in range(N):
            bd, bj = 1e30, 0
            for b in range(N):
                if not used[b] and abs(prev[a]-nxt[b]) < bd:
                    bd, bj = abs(prev[a]-nxt[b]), b
            used[bj] = True; out[a] = nxt[bj]
        return out
    cur = list(np.roots(hG(*pts[0]))); init = list(cur)
    for i in range(1, len(pts)):
        cur = greedy(cur, list(np.roots(hG(*pts[i]))))
    used = [False]*N; perm = [None]*N
    for a in range(N):
        bd, bj = 1e30, -99
        for b in range(N):
            if used[b]: continue
            d = abs(cur[a]-init[b])
            if d < bd: bd, bj = d, b
        used[bj] = True; perm[a] = bj
    assert sorted(perm) == list(range(N))
    return perm, md
perms = []
for name, path in loops.items():
    pp, md = track(path, 2500)
    print(f"  {name:16s} min|D|={md:.2g} perm={pp}")
    perms.append(tuple(pp))
def compose(p, q): return tuple(p[q[i]] for i in range(N))
G = {tuple(range(N))}; stack = list(perms)
while stack:
    g_ = stack.pop()
    for a in perms:
        for hh in (compose(g_, a), compose(a, g_)):
            if hh not in G: G.add(hh); stack.append(hh)
print(f"  |G| = {len(G)} (S5 = 120)")

print("\n=== (C) FORCED SWALLOWTAIL ===")
c2, c3, c4, c5 = symbols("c2 c3 c4 c5")
p1_ = symbols("p1")
pS = p1_*w + c2*w**2 + c3*w**3 + c4*w**4 + c5*w**5
PhiS = integrate(pS, w)
t0 = R(1,2)
conds = [sp.Eq(diff(pS, w).subs(w, t0), 0), sp.Eq(diff(pS, w, 2).subs(w, t0), 0),
         sp.Eq(pS.subs(w, 1), -1), sp.Eq(PhiS.subs(w, 1), 0)]
sol = sp.solve(conds + [sp.Eq(c5, 1)], [p1_, c2, c3, c4, c5], dict=True)[0]
print("  solution (c5 = 1 fixed):", sol)
pSW = sp.expand(pS.subs(sol)); PhiSW = sp.expand(PhiS.subs(sol))
print("  p_swallowtail =", pSW)
print("  checks: p'(1/2) =", diff(pSW, w).subs(w, t0), " p''(1/2) =", diff(pSW, w, 2).subs(w, t0),
      " p'''(1/2) =", diff(pSW, w, 3).subs(w, t0), " p(1) =", pSW.subs(w, 1))
s0 = pSW.subs(w, t0); r0 = t0*s0 - PhiSW.subs(w, t0)
print(f"  swallowtail point: (s,r) = ({s0}, {r0}) = ({float(s0):.6f}, {float(r0):.6f})")
hs = PhiSW - s*s0 + r0
mult = sp.factor(expand(hs - s*0))
print("  h(w; s0, r0) factors:", sp.factor(expand(PhiSW.subs(w, w) - s0*w + r0)))
# escape slope at swallowtail
def pSWn(wv): return float(pSW.subs(w, wv)) if not hasattr(wv, 'imag') else complex(pSW.subs(w, wv))
gam, delt = [], []
for k in np.arange(4, 8.001, 0.5):
    d = 10.0**(-k)
    sv, rv = float(s0) + d, float(r0) + 0.7*d
    roots = np.roots([complex(c) for c in sp.Poly(expand(PhiSW - s*sv + rv), w).all_coeffs()])
    g = [abs(sv - complex(pSW.subs(w, wv))) for wv in roots]
    gam.append(min(g)); delt.append(d*np.hypot(1.0, 0.7))
gam = np.array(gam); delt = np.array(delt)
fit = np.polyfit(np.log(delt), np.log(gam), 1)
print(f"  escape slope at swallowtail: {fit[0]:.4f}  const {np.exp(fit[1]):.4f}   [predict 3/4]")
json.dump({"delt": list(delt), "gam": list(gam), "slope": float(fit[0])},
          open("swallowtail_escape.json","w"))
json.dump({"p": str(pSW), "s0": float(s0), "r0": float(r0)}, open("swallowtail_seed.json","w"))
