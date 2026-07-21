"""
Note 8, stage B: monodromy of F6's 7-sheeted fiber cover. Expect S7 (|G| = 5040),
fold transposition, cusp 3-cycle, acnodes double transpositions, s-loop 6-cycle,
r-loop 7-cycle (w ~ (7s)^{1/6} and (7r)^{1/7} clusters).
"""
import numpy as np, sympy as sp

# Phi6 = -w^7/7 + w^6/6 - 20 w^3/21 + 13 w^2/14 ; h = Phi6 - s w + r
def h_coeffs(s, r):
    return [-1/7, 1/6, 0, 0, -20/21, 13/14, -s, r]

w_, s_, r_ = sp.symbols("w s r")
p6 = -w_**6 + w_**5 - sp.Rational(20,7)*w_**2 + sp.Rational(13,7)*w_
Phi6 = sp.expand(sp.integrate(p6, w_))
D6 = sp.sympify(open("atlas6_wall.txt").read())
disc = sp.lambdify((s_, r_), D6, "numpy")

s_f = float(p6.subs(w_, sp.Rational(1,2)))
r_f = float(sp.Rational(1,2)*p6.subs(w_, sp.Rational(1,2)) - Phi6.subs(w_, sp.Rational(1,2)))
print(f"fold base point: ({s_f}, {r_f})")

CUSP = (0.3043413, 0.0333825)     # refine below
ACN1 = (-1.76416996, 2.29113957)
ACN2 = (-0.10068230, -0.86419633)
# refine the real cusp
wv_c = [complex(v) for v in sp.nroots(sp.diff(p6, w_), n=50)]
for wv in wv_c:
    if abs(wv.imag) < 1e-30:
        sv = complex(p6.subs(w_, wv)); rv = complex(wv*sv - Phi6.subs(w_, wv))
        CUSP = (sv.real, rv.real)
print(f"refined real cusp: {CUSP}")

N = 7
def track(path, steps):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    md = min(abs(disc(s, r)) for s, r in pts)
    def greedy(prev, nxt):
        used = [False]*N; out = [None]*N
        for a in range(N):
            bd, bj = 1e30, 0
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
    assert -99 not in perm and sorted(perm) == list(range(N)), f"POISONED PERM {perm}"
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
    f"fold @({s_f:.5f},{r_f:.5f})": lambda t: (s_f + 0.01*np.exp(1j*t)/np.sqrt(2),
                                                r_f + 1j*0.01*np.exp(1j*t)/np.sqrt(2)),
    "cusp @(0.30434,0.03338)":  lambda t: (CUSP[0] + 0.006*np.exp(1j*t)/np.sqrt(2),
                                            CUSP[1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)),
    "acnode1 @(-1.7642,2.2911)": lambda t: (ACN1[0] + 0.02*np.exp(1j*t)/np.sqrt(2),
                                             ACN1[1] + 1j*0.02*np.exp(1j*t)/np.sqrt(2)),
    "acnode2 @(-0.1007,-0.8642)": lambda t: (ACN2[0] + 0.02*np.exp(1j*t)/np.sqrt(2),
                                              ACN2[1] + 1j*0.02*np.exp(1j*t)/np.sqrt(2)),
    "s=200 e^it, r=1":  lambda t: (200*np.exp(1j*t), 1.0),
    "r=200 e^it, s=0":  lambda t: (0.0, 200*np.exp(1j*t)),
}
perms = {}
for name, path in loops.items():
    p1, md = track(path, 3500)
    p2, _ = track(path, 7000)
    ok = (p1 == p2)
    print(f"{name:34s} min|D6|={md:.3g}  perm={perm_str(p1):24s} refine-agree={ok}", flush=True)
    assert ok, f"loop {name} fails refinement check"
    perms[name] = tuple(p1)

def compose(p, q): return tuple(p[q[i]] for i in range(N))
gens = list(perms.values())
G = {tuple(range(N))}
stack = list(gens)
while stack:
    g = stack.pop()
    if len(G) > 5040:
        raise RuntimeError("GROUP OVERFLOW > 5040 = |S7| -- poison!")
    for a in gens:
        for hh in (compose(g, a), compose(a, g)):
            if hh not in G:
                G.add(hh); stack.append(hh)
print(f"\ngroup generated: |G| = {len(G)} (S7 has 5040)")
ntr = sum(1 for g in G if sum(1 for i in range(N) if g[i] != i) == 2)
print(f"transpositions in G: {ntr} of 21;  |G| == 5040: {len(G) == 5040}")
