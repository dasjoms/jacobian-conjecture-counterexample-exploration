"""
Note 11, stage C: monodromy of F9's 10-sheeted fiber cover -> S10 via JORDAN CERTIFICATE.
Full closure (3.6M) replaced by the rigorous theorem path:
    primitive (transitive on pairs) + contains a transposition  =>  S_n     (C. Jordan)
Machine-checked: (i) fold loop yields a transposition (refinement-verified);
(ii) group transitive; (iii) edge-orbit of a pair under the verified generators
     reaches all 45 pairs -> 2-homogeneous -> primitive.
Plus the usual loop certificates: fold(transposition), cusp(3-cycle),
s-loop(9-cycle w ~ (10s)^{1/9}), r-loop(10-cycle w ~ (10r)^{1/10}).
"""
import numpy as np, sympy as sp

N = 10
def h_coeffs(s, r):
    return [-1/10, 1/9, 0, 0, 0, 0, 0, -44/45, 29/30, -s, r]

w_, s_, r_ = sp.symbols("w s r")
p9 = -w_**9 + w_**8 - sp.Rational(44,15)*w_**2 + sp.Rational(29,15)*w_
Phi10 = sp.expand(sp.integrate(p9, w_))
D9 = sp.sympify(open("atlas9_wall.txt").read())
disc = sp.lambdify((s_, r_), D9, "numpy")

half = sp.Rational(1,2)
s_f = float(p9.subs(w_, half)); r_f = float(half*p9.subs(w_, half) - Phi10.subs(w_, half))
CUSPS = []
for wv in sp.nroots(sp.diff(p9, w_), n=60):
    if abs(sp.im(wv)) < sp.Float("1e-50"):
        sv = complex(p9.subs(w_, wv)); rv = complex(wv*sv - Phi10.subs(w_, wv))
        CUSPS.append((sv.real, rv.real, float(sp.re(wv))))
print(f"fold ({s_f:.5f},{r_f:.5f});  real cusps: {[(round(a,5), round(b,5)) for a,b,_ in CUSPS]}")
CU = [c for c in CUSPS if c[2] > 0][0]

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
            used[bj] = True; return_out_bj = bj; out[a] = nxt[bj]
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
            dd = abs(cur[a]-init[b])
            if dd < bd: bd, bj = dd, b
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
    f"cusp @({CU[0]:.5f},{CU[1]:.5f})": lambda t: (CU[0] + 0.006*np.exp(1j*t)/np.sqrt(2),
                                                    CU[1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)),
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
    print(f"{name:38s} min|D9|={md:.3g}  perm={perm_str(p1):36s} refine-agree={ok}", flush=True)
    assert ok, f"refinement mismatch on {name}"
    perms[name] = tuple(p1)

# ---- Jordan certificate ----
gens = list(perms.values())
def compose(p, q): return tuple(p[q[i]] for i in range(N))
# transitivity (points)
orb = {0}; stack = [0]
while stack:
    v = stack.pop()
    for g in gens:
        nv = g[v]
        if nv not in orb: orb.add(nv); stack.append(nv)
transitive = (len(orb) == N)
# 2-homogeneity: orbit of the pair {0,1}
edge = frozenset((0,1)); seen = {edge}; stack = [edge]
while stack:
    e = stack.pop()
    a, b = tuple(e)
    for g in gens:
        ne = frozenset((g[a], g[b]))
        if ne not in seen: seen.add(ne); stack.append(ne)
two_homog = (len(seen) == N*(N-1)//2)
# transposition present
tr = [g for g in gens if sum(1 for i in range(N) if g[i] != i) == 2]
print(f"\nJORDAN CERTIFICATE:")
print(f"  transitive: {transitive}   2-homogeneous (all {N*(N-1)//2} pairs reached): {two_homog}")
print(f"  transposition among verified loops: {[perm_str(t) for t in tr]}")
verdict = transitive and two_homog and len(tr) > 0
print(f"  => primitive + transposition => G = S{N} (Jordan).   VERDICT: {verdict}")
open("atlas9_monodromy.txt","w").write(
    f"jordan_transitive={transitive} two_homogeneous={two_homog} transposition={[perm_str(t) for t in tr]}\n" +
    "\n".join(f"{k}: {perm_str(v)}" for k, v in perms.items()))

# bonus: partial closure (cap 400k) to exhibit a lower bound on |G|
import time
t0 = time.time()
G = set(gens) | {tuple(range(N))}; stack = list(gens)
while stack and len(G) < 400000:
    g = stack.pop()
    for a in gens:
        for hh in (compose(g, a), compose(a, g)):
            if hh not in G:
                G.add(hh); stack.append(hh)
ntr = sum(1 for g in G if sum(1 for i in range(N) if g[i] != i) == 2)
print(f"bonus: partial closure |G| >= {len(G)} [{time.time()-t0:.0f}s], transpositions so far: {ntr}")
