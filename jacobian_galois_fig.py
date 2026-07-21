#!/usr/bin/env python3
"""NOTE 19 figure v2: THE PIN'S TRANSPOSITION (4 panels, fixed geometry)."""
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mpmath as mp

w = sp.symbols('w')
def tower(d):
    c0 = sp.Rational(6, d * (d + 1))
    p = sp.expand(2*w - 3*w**2 + w*(1-w)*(w**(d-2) - c0))
    Phi = sp.expand(sp.integrate(p, w))
    return p, Phi

fig = plt.figure(figsize=(13.6, 10.8))
fig.patch.set_facecolor('#fdfcf8')
gs = fig.add_gridspec(2, 2, left=0.07, right=0.97, top=0.90, bottom=0.08, wspace=0.24, hspace=0.34)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1]); axB.axis('off')
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])
fig.suptitle('NOTE 19 — THE PIN’S TRANSPOSITION: exact Galois certificates for the atlas',
             fontsize=13.5)

# ---------- (a) two loops around the pin swap exactly two sheets ----------
d = 4
p_, Phi_ = tower(d)
pol = sp.Poly(Phi_, w)
n = d + 1
base = [mp.mpf(str(c)) for c in pol.all_coeffs()]          # high -> low
thetas = np.linspace(0, 2*np.pi, 420, endpoint=False)
sel = []; prev = None
mp.mp.dps = 40
for th in thetas:
    s0 = -1 + 0.05*np.exp(1j*th)
    coefs = list(base)
    coefs[n - 1] = coefs[n - 1] - s0
    coefs[n] = coefs[n] - 1
    rt = list(mp.polyroots(coefs, maxsteps=250))
    pair = sorted(rt, key=lambda z: abs(z - 1))[:2]
    if prev is not None:
        cands = [(pair[0], pair[1]), (pair[1], pair[0])]
        pair = list(min(cands, key=lambda cand: abs(cand[0]-prev[0]) + abs(cand[1]-prev[1])))
    prev = pair
    sel.append(pair)
A = np.array([s_[0] for s_ in sel], dtype=complex)
B = np.array([s_[1] for s_ in sel], dtype=complex)
axA.plot(A.real, A.imag, '.', color='#d62728', ms=3.2)
axA.plot(B.real, B.imag, '.', color='#1f77b4', ms=3.2)
axA.plot([1], [0], 'k*', ms=18)
axA.plot([], [], color='#d62728', lw=2, label='sheet A')
axA.plot([], [], color='#1f77b4', lw=2, label='sheet B')
axA.set_xlim(0.3, 2.05); axA.set_ylim(-1.15, 1.15)
axA.set_aspect('equal'); axA.grid(alpha=0.25)
axA.text(0.36, 0.92, 'one loop around the pin:\nA and B exchange partners;\nsecond loop returns them\n⟹ local monodromy is a\nsingle transposition', fontsize=9,
         bbox=dict(fc='white', ec='#999', boxstyle='round,pad=0.3'))
axA.set_xlabel('Re w'); axA.set_ylabel('Im w')
axA.set_title('(a) C1 at d=4: loop  s = −1 + 0.05·e^{iθ}, r = −1\naround the pin; sheets tracked at 40 digits', fontsize=10.5)
axA.legend(fontsize=8.5, loc='upper right')

# ---------- (b) Newton-polygon schematics ----------
axB.set_xlim(0, 10); axB.set_ylim(0, 10)
axB.set_title('(b) C2 + C3: inertia by exact Newton–Puiseux\nslopes (d = 5, n = 6; verified symbolically, d = 2..12)', fontsize=10.5)
axB.add_patch(Rectangle((0.35, 4.2), 4.6, 5.3, fc='#eaf3ff', ec='#225588', lw=1.5))
axB.text(2.65, 9.05, 'at s = ∞: edge slope 1/(n−1)', ha='center', fontsize=9.5, color='#17456e')
pts = [(0, 0), (1, -1)] + [(i, 0) for i in range(2, 7)]
xs = [1.0 + q*0.5 for q, _ in pts]; ys = [8.35 + v*2.9 for _, v in pts]
axB.plot(xs, ys, 'o', color='#17456e', ms=6)
axB.plot([xs[0], xs[1]], [ys[0], ys[1]], '-', color='#999', lw=1)
axB.plot([xs[1], xs[-1]], [ys[1], ys[-1]], '-', color='#d62728', lw=2.4)
axB.text(2.62, 4.55, 'G(W,0) = aₙWⁿ − W,\nWⁿ⁻¹ = 1/aₙ:\nn−1 single-orbit roots\n⟹ inertia ⊇ (n−1)-cycle', fontsize=8.8, ha='center',
         bbox=dict(fc='white', ec='#d62728', boxstyle='round,pad=0.3'))
axB.add_patch(Rectangle((5.2, 4.2), 4.6, 5.3, fc='#fff0ec', ec='#bb3311', lw=1.5))
axB.text(7.5, 9.05, 'at r = ∞: edge slope 1/n', ha='center', fontsize=9.5, color='#8c2208')
pts = [(0, -1)] + [(i, 0) for i in range(1, 7)]
xs = [5.85 + q*0.5 for q, _ in pts]; ys = [8.35 + v*2.9 for _, v in pts]
axB.plot(xs, ys, 'o', color='#8c2208', ms=6)
axB.plot([xs[0], xs[-1]], [ys[0], ys[-1]], '-', color='#bb3311', lw=2.4)
axB.text(7.5, 4.55, 'G(W,0) = aₙWⁿ + 1,\nWⁿ = −1/aₙ:\nn single-orbit roots\n⟹ inertia ⊇ n-cycle', fontsize=8.8, ha='center',
         bbox=dict(fc='white', ec='#bb3311', boxstyle='round,pad=0.3'))
axB.text(5.0, 2.95, 'the middle coefficients sit strictly above each edge (weight-0\nvaluations) — symbolic, and formula-true for every d', ha='center', fontsize=8.8, color='#444')
axB.text(5.0, 1.15, 'group lemma: n-cycle ⟹ transitive; (n−1)-cycle (fixed point moved\naround by conjugates of the n-cycle) ⟹ every point stabilizer transitive ⟹ 2-transitive;\n2-transitive + one transposition ⟹ ALL transpositions ⟹ Sₙ', ha='center', fontsize=8.6,
         color='#1c6b1c', bbox=dict(fc='#f4fbf4', ec='#2ca02c', boxstyle='round,pad=0.35'))

# ---------- (c) d = 3 coalescence ----------
p3, Phi3 = tower(3)
tt = np.linspace(-3.6, 2.6, 6000)
sv = np.array([float(p3.subs(w, t_)) for t_ in tt])
rv = np.array([float(t_*p3.subs(w, t_) - Phi3.subs(w, t_)) for t_ in tt])
axC.plot(sv, rv, color='#225588', lw=1.3, alpha=0.85, label='wall of Φ₃, parametrized by t')
axC.plot([-1], [-1], 'ko', ms=9)
axC.annotate('the pin (−1,−1)\n t = 1 maps here', xy=(-1.02, -1.06), xytext=(-3.4, -3.0), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='k'))
axC.plot([-2.5], [float(-2.5*p3.subs(w,-2.5)-Phi3.subs(w,-2.5))], 'v', ms=10, mfc='none', mec='#8c2208')
axC.annotate('t = −2 maps here TOO\ncorner–shadow coalescence (note 17)\npin gcd = (w−1)(w+2)\n⟹ C1 moves to t₀ = 2',
             xy=(-0.99, -0.9), xytext=(-2.9, 0.6), fontsize=9,
             arrowprops=dict(arrowstyle='->', color='#8c2208'))
s2 = float(p3.subs(w, 2)); r2 = float(2*p3.subs(w, 2) - Phi3.subs(w, 2))
axC.plot([s2], [r2], 's', ms=10, color='#1c6b1c')
axC.annotate('t₀ = 2: clean anchor\n(single double sheet)', xy=(s2, r2), xytext=(-0.4, 1.2), fontsize=9, color='#1c6b1c',
             arrowprops=dict(arrowstyle='->', color='#1c6b1c'))
axC.set_xlim(-3.8, 2.9); axC.set_ylim(-4.2, 2.1)
axC.grid(alpha=0.25); axC.legend(fontsize=8.5, loc='lower left')
axC.set_xlabel('s'); axC.set_ylabel('r')
axC.set_title('(c) d = 3 honesty: the pin IS a 2-branch point there —\nthe certificate reroutes to another smooth wall point', fontsize=10.5)

# ---------- (d) coverage matrix ----------
rows = list(range(2, 13))
cols = ['C1\ntransp.', 'C2\n(n−1)-cyc', 'C3\nn-cyc', 'Gal=Sₙ\ncertified', 'Dedekind\ntrinity']
grid = np.full((len(rows), len(cols)), 2.0)
for i, d in enumerate(rows):
    if d > 7:
        grid[i, 4] = 1.0
    if d == 3:
        grid[i, 0] = 1.75
cmap = matplotlib.colors.ListedColormap(['#dddddd', '#ffe9b8', '#b8e3b8'])
axD.imshow(grid, cmap=cmap, vmin=0, vmax=2, aspect='auto')
axD.set_xticks(range(len(cols))); axD.set_xticklabels(cols, fontsize=8.6)
axD.set_yticks(range(len(rows))); axD.set_yticklabels([f'd={d} (S{d+1})' for d in rows], fontsize=8.6)
for i, d in enumerate(rows):
    for j in range(len(cols)):
        if j == 4 and d > 7:
            txt, col = 'n, n−1\nonly', '#8a6d1a'
        elif j == 0 and d == 3:
            txt, col = '✓ (t₀=2)', '#1c6b1c'
        else:
            txt, col = '✓', '#1c6b1c'
        axD.text(j, i, txt, ha='center', va='center', fontsize=8.2, color=col)
axD.set_xlim(-0.5, 4.5); axD.set_ylim(len(rows)-0.5, -0.65)
axD.set_title('(d) every atlas chamber: braid now PROVEN full,\nplus all three exact certificates hold to d = 30', fontsize=10.5)
axD.text(2, len(rows)-0.15, 'dedekind d≥8: transposition-type density 1/(2(n−2)!) ≤ 1e-4\n— corroboration only; the proof never needs it',
         ha='center', va='top', fontsize=7.8, color='#555')

fig.savefig('/home/user/galois_figure.png', dpi=170)
print('saved galois_figure.png')
