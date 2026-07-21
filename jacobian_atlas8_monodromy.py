"""
Note 10, stage C: monodromy of F8's 9-sheeted fiber cover. Expect S9 (|G| = 362880).
Loops: fold @w=1/2 (transposition), real cusp (3-cycle), acnode, s-loop (8-cycle), r-loop (9-cycle).
Poison guards: disjoint-support -1 perm detector (assert sorted(perm)==range(9)), overflow cap.
"""
import numpy as np, sympy as sp
from collections import Counter

N = 9
def h_coeffs(s, r):
    return [-1/9, 1/8, 0, 0, 0, 0, -35/36, 23/24, -s, r]

w_, s_, r_ = sp.symbols("w s r")
p8 = -w_**8 + w_**7 - sp.Rational(35,12)*w_**2 + sp.Rational(23,12)*w_
Phi9 = sp.expand(sp.integrate(p8, w_))
D8 = sp.sympify(open("atlas8_wall.txt").read())
disc = sp.lambdify((s_, r_), D8, "numpy")

half = sp.Rational(1,2)
s_f = float(p8.subs(w_, half)); r_f = float(half*p8.subs(w_, half) - Phi9.subs(w_, half))
CUSP = None
for wv in sp.nroots(sp.diff(p8, w_), n=60):
    if abs(sp.im(wv)) < sp.Float("1e-50"):
        sv = complex(p8.subs(w_, wv)); rv = complex(wv*sv - Phi9.subs(w_, wv))
        CUSP = (sv.real, rv.real); print(f"real cusp contact t={complex(wv):.6f}")
print(f"fold ({s_f:.5f},{r_f:.5f})  real cusp {CUSP}")
ACN = None
dd = sp.sympify(open("atlas8_bitangents.json").read()) if False else None

def track(path, steps):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    md = min(abs(disc(s, r)) for s, r in pts)
    def greedy(prev, nxt):
        used = [False]*N; out = [None]*N
        for a in range(N):
            bd, bj = 1e30, -1
            for b in range(N):
                if not used[b] and abs(prev[a]-nxt[b]) < bd:
                    bd, bj = abs(prev[a]-nxt[b]), b
            used[bj] = True; out[a] = nxt[bj]
        return out
    cur = list(np.roots(h_coeffs(*pts[0])))
    init = list(cur)
    for i in range(1, len(pts)):
        cur = greedy(cur, list(np.roots(h_coeffs(*pts[i]))))
    used = [False]*N; perm = [None]*N
    for a in range(N):
        bd, bj = 1e30, -99
        for b in range(N):
            if used[b]: continue
            d = abs(cur[a]-init[b])
            if d < bd: bd, bj = d, b
        used[bj] = True; perm[a] = bj
    assert sorted(perm) == list(range(N)), f"POISONED {perm}"
    return perm, md

def perm_str(pp):
    seen = [False]*N; cyc = []
    for i in range(N):
        if seen[i]: continue
        j, c = i, []
        while not seen[j]:
            seen[j] = True; c.append(j+1); j = pp[j]
        if len(c) > 1: cyc.append(tuple(c))
    return " x ".join("("+",".join(map(str, c))+")" for c in cyc) or "id"

loops = {
    f"fold @({s_f:.5f},{r_f:.5f})": lambda t: (s_f + 0.008*np.exp(1j*t)/np.sqrt(2),
                                                r_f + 1j*0.008*np.exp(1j*t)/np.sqrt(2)),
    f"cusp @({CUSP[0]:.5f},{CUSP[1]:.5f})": lambda t: (CUSP[0] + 0.006*np.exp(1j*t)/np.sqrt(2),
                                                        CUSP[1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)),
    "s=200 e^it, r=1":  lambda t: (200*np.exp(1j*t), 1.0),
    "r=200 e^it, s=0":  lambda t: (0.0, 200*np.exp(1j*t)),
}
perms = {}
for name, path in loops.items():
    p1, md = track(path, 4000)
    p2, _ = track(path, 8000)
    ok = (p1 == p2)
    if not ok:
        p3, _ = track(path, 16000)
        ok = (p2 == p3); p1 = p2
    print(f"{name:38s} min|D8|={md:.3g}  perm={perm_str(p1):30s} refine-agree={ok}", flush=True)
    assert ok, f"refinement mismatch on {name}"
    perms[name] = tuple(p1)

def compose(p, q): return tuple(p[q[i]] for i in range(N))
gens = list(perms.values())
G = set(gens) | {tuple(range(N))}
stack = list(gens)
import time
t0 = time.time()
while stack:
    g = stack.pop()
    if len(G) > 362880:
        raise RuntimeError("GROUP OVERFLOW > 9! = 362880")
    for a in gens:
        for hh in (compose(g, a), compose(a, g)):
            if hh not in G:
                G.add(hh); stack.append(hh)
print(f"\ngroup generated: |G| = {len(G)}  (S9 has 362880)   [{time.time()-t0:.0f}s]")
ntr = sum(1 for g in G if sum(1 for i in range(N) if g[i] != i) == 2)
print(f"transpositions in G: {ntr} of 36;  |G| == 362880: {len(G) == 362880}")
open("atlas8_monodromy.txt","w").write(
    f"|G|={len(G)} transpositions={ntr}\n" +
    "\n".join(f"{k}: {perm_str(v)}" for k, v in perms.items()))
