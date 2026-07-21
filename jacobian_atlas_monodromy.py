"""
Stage B (fixed): monodromy of F4's 5-sheeted fiber family in COMPLEX (s,r)-space.
Each loop lies in a complex line through the wall point / around infinity and is
verified to keep |Disc| bounded below (no wall crossing).
Roots: isolated per step via np.roots and greedy nearest-neighbor continuation.
Fine-step consistency check: every loop recomputed at 2x steps must agree.
"""
import numpy as np

C5, C4, C3, C2 = -1/5, 1/4, -9/10, 17/20
def h_coeffs(s, r):
    return [C5, C4, C3, C2, -s, r]

D4_terms = [(-20000000, 4, 0), (20000000, 3, 1), (18900000, 3, 0), (-57100000, 2, 2),
            (75376000, 2, 1), (-46613461, 2, 0), (38210000, 1, 3), (-67741050, 1, 2),
            (45893931, 1, 1), (-6657115, 1, 0), (-8192000, 0, 5), (17058675, 0, 4),
            (-12278715, 0, 3), (1957975, 0, 2)]
def disc_val(s, r):
    return sum(c * r**rr * s**ss for c, rr, ss in D4_terms)

def track(path, steps, check_wall=True):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    if check_wall:
        md = min(abs(disc_val(s, r)) for s, r in pts)
    else:
        md = None
    def greedy(prev, nxt):
        used = [False]*len(prev); out = [None]*len(prev)
        for a in range(len(prev)):
            bd, bj = 1e30, 0
            for b in range(len(nxt)):
                if not used[b] and abs(prev[a]-nxt[b]) < bd:
                    bd, bj = abs(prev[a]-nxt[b]), b
            used[bj] = True; out[a] = nxt[bj]
        return out
    cur = list(np.roots(h_coeffs(*pts[0])))
    init = list(cur)
    for i in range(1, len(pts)):
        cur = greedy(cur, list(np.roots(h_coeffs(*pts[i]))))
    # final permutation: match cur (at theta=2pi-eps, near pts[0]) to init
    used = [False]*5; perm = [None]*5
    for a in range(5):
        bd, bj = 1e30, 0
        for b in range(5):
            if used[b]: continue
            d = abs(cur[a]-init[b])
            if d < bd: bd, bj = d, b
        used[bj] = True; perm[a] = bj
    return perm, md

def perm_str(pp):
    seen = [False]*5; cyc = []
    for i in range(5):
        if seen[i]: continue
        j, c = i, []
        while not seen[j]:
            seen[j] = True; c.append(j+1); j = pp[j]
        if len(c) > 1: cyc.append(tuple(c))
    return " x ".join("("+",".join(map(str, c))+")" for c in cyc) or "id"

loops = {
    "fold @(-1,-1)":      lambda t: (-1 + 0.03*np.exp(1j*t)/np.sqrt(2), -1 + 1j*0.03*np.exp(1j*t)/np.sqrt(2)),
    "cusp @(0.2921,0.034)": lambda t: (0.2921225227866996 + 0.008*np.exp(1j*t)/np.sqrt(2),
                                       0.0340042595568076877 + 1j*0.008*np.exp(1j*t)/np.sqrt(2)),
    "node @(0.9842,0.6762)": lambda t: (0.9841655956599125 + 0.01*np.exp(1j*t)/np.sqrt(2),
                                        0.67623935 + 1j*0.01*np.exp(1j*t)/np.sqrt(2)),
    "s=50 e^it, r=1":     lambda t: (50*np.exp(1j*t), 1.0),
    "r=50 e^it, s=0":     lambda t: (0.0, 50*np.exp(1j*t)),
}
perms = {}
for name, path in loops.items():
    p1, md = track(path, 2500)
    p2, _ = track(path, 5000)
    ok = (p1 == p2)
    print(f"{name:28s}  min|Disc|={md:.3g}  perm={perm_str(p1):15s} refine-agree={ok}")
    assert ok, f"loop {name} fails refinement check"
    perms[name] = tuple(p1)

def compose(p, q): return tuple(p[q[i]] for i in range(5))
gens = list(perms.values())
G = {tuple(range(5))}
stack = list(gens)
while stack:
    g = stack.pop()
    for a in gens:
        for h in (compose(g, a), compose(a, g)):
            if h not in G:
                G.add(h); stack.append(h)
print(f"\ngroup generated: |G| = {len(G)} (S5 has 120)")
nt = [g for g in G if sum(1 for i in range(5) if g[i] != i) == 2]
print(f"transpositions found: {len(nt)};  |G| == 120: {len(G) == 120}")
