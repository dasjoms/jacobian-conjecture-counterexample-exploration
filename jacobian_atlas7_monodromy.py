"""
Note 9, stage B: monodromy of F7's 8-sheeted fiber cover. Expect S8 (|G| = 40320).
"""
import numpy as np, sympy as sp

def h_coeffs(s, r):
    return [-1/8, 1/7, 0, 0, 0, -27/28, 53/56, -s, r]

w_, s_, r_ = sp.symbols("w s r")
p7 = -w_**7 + w_**6 - sp.Rational(81,28)*w_**2 + sp.Rational(53,28)*w_
Phi7 = sp.expand(sp.integrate(p7, w_))
D7 = sp.sympify(open("atlas7_wall.txt").read())
disc = sp.lambdify((s_, r_), D7, "numpy")

half = sp.Rational(1,2)
s_f = float(p7.subs(w_, half)); r_f = float(half*p7.subs(w_, half) - Phi7.subs(w_, half))
CUSP = None
for wv in sp.nroots(sp.diff(p7, w_), n=50):
    if abs(sp.im(wv)) < sp.Float("1e-40") and 0 < float(sp.re(wv)) < 1:
        sv = complex(p7.subs(w_, wv)); rv = complex(wv*sv - Phi7.subs(w_, wv))
        CUSP = (sv.real, rv.real)
print(f"fold ({s_f:.5f},{r_f:.5f})  cusp {CUSP}")

N = 8
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
    assert -99 not in perm and sorted(perm) == list(range(N)), f"POISONED {perm}"
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

CRU = (-0.90945711, -0.91031417)
ACN = (0.13989048, 1.54077491)
loops = {
    f"fold @({s_f:.5f},{r_f:.5f})": lambda t: (s_f + 0.008*np.exp(1j*t)/np.sqrt(2),
                                                r_f + 1j*0.008*np.exp(1j*t)/np.sqrt(2)),
    "cusp @(0.31048,0.03400)":  lambda t: (CUSP[0] + 0.006*np.exp(1j*t)/np.sqrt(2),
                                            CUSP[1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)),
    "crunode @(-0.9095,-0.9103)": lambda t: (CRU[0] + 0.02*np.exp(1j*t)/np.sqrt(2),
                                              CRU[1] + 1j*0.02*np.exp(1j*t)/np.sqrt(2)),
    "acnode @(0.1399,1.5408)":  lambda t: (ACN[0] + 0.03*np.exp(1j*t)/np.sqrt(2),
                                            ACN[1] + 1j*0.03*np.exp(1j*t)/np.sqrt(2)),
    "s=200 e^it, r=1":  lambda t: (200*np.exp(1j*t), 1.0),
    "r=200 e^it, s=0":  lambda t: (0.0, 200*np.exp(1j*t)),
}
perms = {}
for name, path in loops.items():
    p1, md = track(path, 4000)
    p2, _ = track(path, 8000)
    ok = (p1 == p2)
    print(f"{name:34s} min|D7|={md:.3g}  perm={perm_str(p1):26s} refine-agree={ok}", flush=True)
    assert ok
    perms[name] = tuple(p1)

def compose(p, q): return tuple(p[q[i]] for i in range(N))
gens = list(perms.values())
G = {tuple(range(N))}
stack = list(gens)
while stack:
    g = stack.pop()
    if len(G) > 40320:
        raise RuntimeError("GROUP OVERFLOW > 40320!")
    for a in gens:
        for hh in (compose(g, a), compose(a, g)):
            if hh not in G:
                G.add(hh); stack.append(hh)
print(f"\ngroup generated: |G| = {len(G)} (S8 has 40320)")
ntr = sum(1 for g in G if sum(1 for i in range(N) if g[i] != i) == 2)
print(f"transpositions in G: {ntr} of 28;  |G| == 40320: {len(G) == 40320}")
