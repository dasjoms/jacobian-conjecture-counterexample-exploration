"""
Note 7, stage C second addendum:
1. hinge x-pair matching: full sorted eps-sweep x-list vs exact hinge roots.
2. node preimages at 100-digit mpmath precision: coalesced pairs, 2 bounded + 4 escaping.
3. hinge antipodal-pair theorem: symbolic identity check (q2 = p1/2, involution u -> -2/p1 - u).
4. grid of real-root counts over (s,r) for the region atlas figure.
"""
import numpy as np, json
from mpmath import mp, mpf, mpc, polyroots

mp.dps = 100

def p5w(w):
    return -w**5 + w**4 - 14*w**2/5 + 9*w/5

print("== 1. hinge pairs matched against eps-sweep ==")
def hinge(A, B):
    qa = A*81/25 - B**2*9/10; qb = A*18/5 - B**2; qc = A
    us = np.roots([qa, qb, qc])
    return sorted([complex(-(1+1.8*u)/B) for u in us], key=lambda z: z.real)
for (A, B) in [(2.0, 3.0), (-0.5483, -1.7812)]:
    xs = hinge(A, B)
    roots = np.roots([-1/6, 1/5, 0, -14/15, 9/10, -B*1e-5, A*1e-10])
    gam = [B*1e-5 - p5w(wv) for wv in roots]
    xx = sorted([1e-5/g for g in gam], key=lambda z: abs(z))
    print(f"  (A,B)=({A},{B}): hinge {[(round(v.real,5), round(v.imag,5)) for v in xs]}")
    print(f"    eps-sweep all |x|: {[(round(v.real,5), round(v.imag,5)) for v in xx]}")
    d = min(abs(xs[0]-xx[-1]) + abs(xs[1]-xx[-2]), abs(xs[0]-xx[-2]) + abs(xs[1]-xx[-1]))
    print(f"    hinge-vs-two-largest mismatch: {d:.2e}")

print("\n== 2. node fibers at 100 digits ==")
nodes = json.load(open("atlas5_bitangents.json"))
for a, b, sv, rv in nodes["nodes"]:
    sv, rv = complex(sv), complex(rv)
    if not (abs(sv.imag) < 1e-10 and abs(rv.imag) < 1e-10):
        continue  # real nodes only (crunode + acnode)
    coeffs = [mpf(-1)/6, mpf(1)/5, 0, mpf(-14)/15, mpf(9)/10, mpc(-sv.real), mpc(rv.real)]
    roots = polyroots(coeffs, maxsteps=500)
    gams = [mpf(sv.real) - p5w(mpc(z)) for z in roots]
    esc = sum(1 for g in gams if abs(g) < mpf(10)**(-40))
    print(f"  node ({sv.real:.8f},{rv.real:.8f}): |gamma| values sorted:", flush=True)
    for g in sorted(gams, key=lambda gg: abs(gg)):
        print(f"    |gamma| = {abs(g):.3e}")
    print(f"    => escaping sheets {esc}, bounded preimages {6-esc}")

print("\n== 3. hinge antipodal theorem (symbolic) ==")
import sympy as sp
u, p1 = sp.symbols("u p1")
q2 = p1/2
uh = -2/p1 - u
lhs = u*(1+q2*u); rhs = sp.expand(uh*(1+q2*uh))
print("  u(1+q2 u) invariant under involution:", sp.simplify(lhs - rhs) == 0)
print("  1+p1 u flips sign:", sp.simplify(1 + p1*uh + (1 + p1*u)) == 0)

print("\n== 4. real-root-count grid (200x200 over [-4,4]^2) ==")
n = 200
S = np.linspace(-4, 4, n); R = np.linspace(-4, 4, n)
grid = np.zeros((n, n), dtype=int)
for i, r in enumerate(R):
    for j, sv in enumerate(S):
        roots = np.roots([-1/6, 1/5, 0, -14/15, 9/10, -sv, r])
        grid[i, j] = sum(1 for wv in roots if abs(wv.imag) < 1e-7)
vals, counts = np.unique(grid, return_counts=True)
print("  grid counts:", dict(zip(vals.tolist(), counts.tolist())))
np.savez("atlas5_grid.npz", S=S, R=R, grid=grid)
