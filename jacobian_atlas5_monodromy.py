"""
Note 7, stage B: monodromy of F5's 6-sheeted fiber cover in complex (s,r)-space.
Same protocol as note 6: complex-line loops, min|D5| no-cross certificate,
2x-refinement agreement assertion, closure with permutation-sanity guard.
Predictions: fold -> transposition; cusp -> 3-cycle; nodes -> double transpositions;
s-loop -> 5-cycle (w ~ (6s)^{1/5}); r-loop -> 6-cycle (w ~ (6r)^{1/6}); group S6, |G|=720.
"""
import numpy as np

# Phi5 = -w^6/6 + w^5/5 - 14 w^3/15 + 9 w^2/10 ; h = Phi5 - s w + r
def h_coeffs(s, r):
    return [-1/6, 1/5, 0, -14/15, 9/10, -s, r]

D5_terms = [(3037500000,5,0), (-3037500000,4,1), (11907000000,4,0), (-50625000,3,2),
            (-19265850000,3,1), (15615855000,3,0), (9797625000,2,3), (8794980000,2,2),
            (-23436459000,2,1), (8496467856,2,0), (-6968981250,1,4), (2780676000,1,3),
            (8071107300,1,2), (-6147920736,1,1), (892142910,1,0), (1220703125,0,6),
            (-46159500,0,5), (-2098753050,0,4), (1564980260,0,3), (-247817475,0,2)]
def disc_val(s, r):
    return sum(c * r**rr * s**ss for c, rr, ss in D5_terms)

N = 6
def track(path, steps):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    md = min(abs(disc_val(s, r)) for s, r in pts)
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

# fold point from parametrization at t = 1/2:
s_f = 0.23125; r_f = 0.5*s_f - (-(0.5)**6/6 + (0.5)**5/5 - 14*(0.5)**3/15 + 9*(0.5)**2/10)
loops = {
    "fold @(0.23125,%0.5f)" % r_f: lambda t: (s_f + 0.01*np.exp(1j*t)/np.sqrt(2),
                                               r_f + 1j*0.01*np.exp(1j*t)/np.sqrt(2)),
    "cusp @(0.29716,0.03303)":  lambda t: (0.2971583157757676 + 0.008*np.exp(1j*t)/np.sqrt(2),
                                            0.03302630284803969 + 1j*0.008*np.exp(1j*t)/np.sqrt(2)),
    "crunode @(-0.8819,-0.8834)": lambda t: (-0.88188341 + 0.01*np.exp(1j*t)/np.sqrt(2),
                                             -0.88337174 + 1j*0.01*np.exp(1j*t)/np.sqrt(2)),
    "acnode @(0.4920,-0.4641)": lambda t: (0.49202797 + 0.01*np.exp(1j*t)/np.sqrt(2),
                                           -0.46406377 + 1j*0.01*np.exp(1j*t)/np.sqrt(2)),
    "s=50 e^it, r=1":    lambda t: (50*np.exp(1j*t), 1.0),
    "r=200 e^it, s=0":   lambda t: (0.0, 200*np.exp(1j*t)),
}
perms = {}
for name, path in loops.items():
    p1, md = track(path, 3000)
    p2, _ = track(path, 6000)
    ok = (p1 == p2)
    print(f"{name:32s} min|D5|={md:.3g}  perm={perm_str(p1):22s} refine-agree={ok}", flush=True)
    assert ok, f"loop {name} fails refinement check"
    perms[name] = tuple(p1)

def compose(p, q): return tuple(p[q[i]] for i in range(N))
gens = list(perms.values())
G = {tuple(range(N))}
stack = list(gens)
while stack:
    g = stack.pop()
    if len(G) > 720:
        raise RuntimeError("GROUP OVERFLOW > 720 = |S6| -- poisoned generator slipped through!")
    for a in gens:
        for hh in (compose(g, a), compose(a, g)):
            if hh not in G:
                G.add(hh); stack.append(hh)
print(f"\ngroup generated: |G| = {len(G)} (S6 has 720)")
ntr = sum(1 for g in G if sum(1 for i in range(N) if g[i] != i) == 2)
print(f"transpositions in G: {ntr} of {N*(N-1)//2};  |G| == 720: {len(G) == 720}")
