#!/usr/bin/env python3
"""NOTE 18 figure: THE MISSING BRAID (4 panels)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle

def Epoly(v):
    return -27*v**3 - 94.5*v**2 - 111*v - 44.5

fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.6))
fig.patch.set_facecolor('#fdfcf8')

# ---------------- (a) the canal: log-log slope -3 ----------------
ax = axes[0, 0]; ax.set_facecolor('#fdfcf8')
xs = np.logspace(-1.3, 0.9, 260)
slopes = []
for v0, col in [(0.5, '#1f77b4'), (1.0, '#d62728'), (1.7, '#2ca02c'), (2.6, '#9467bd')]:
    ys = v0 / xs
    det = np.abs(ys**3 * Epoly(xs * ys))
    ax.loglog(xs, det, color=col, lw=2.2, label=f'v = xy = {v0}')
    p = np.polyfit(np.log(xs), np.log(det), 1)
    slopes.append(p[0])
xm = np.logspace(-1.3, 0.9, 50)
ax.loglog(xm, 8.2 * xm**-3, 'k--', lw=1.1, alpha=0.55)
ax.text(0.105, 8.2 * 0.105**-3 * 1.25, 'x$^{-3}$ guide', fontsize=9, style='italic', color='k', alpha=0.7)
ax.set_xlabel('|x| (moving along a hyperbola xy = v)', fontsize=10.5)
ax.set_ylabel('|det JF$_2$|  =  |y$^3$ E(v)|', fontsize=10.5)
ax.text(0.03, 0.955, f'fitted slopes: {slopes[0]:.9f} (×4: all equal)', transform=ax.transAxes,
        fontsize=9, va='top', bbox=dict(fc='white', ec='#888', boxstyle='round,pad=0.3'))
ax.set_title('(a) THE CANAL:  det = y³ E(xy)\nexact slope −3 on every hyperbola xy = v: never constant', fontsize=10.5)
ax.legend(fontsize=9, loc='lower left'); ax.grid(alpha=0.25)

# ---------------- (b) braid ledger (schematic) ----------------
ax = axes[0, 1]; ax.axis('off')
ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title('(b) THE BRAID LEDGER:\nweight-0 monomial rings, machine-checked', fontsize=10.5)
# left: 3D
ax.add_patch(Rectangle((0.45, 5.6), 4.1, 3.9, fc='#eaf3ff', ec='#225588', lw=1.6))
ax.text(2.5, 9.05, 'C³, weights (1,−1,−2)', ha='center', fontsize=10.5, fontweight='bold', color='#17456e')
ax.text(2.5, 8.35, 'weight-0 monomials:\n{A−B−2C=0} = {v$^B$ t$^C$}', ha='center', fontsize=9.5,
        va='top', bbox=dict(fc='white', ec='#225588', boxstyle='round,pad=0.35'))
ax.text(1.25, 6.15, 'v = xy', fontsize=11, ha='center', color='#17456e')
ax.text(3.75, 6.15, 't = x²z', fontsize=11, ha='center', color='#17456e')
ax.annotate('', xy=(2.5, 4.9), xytext=(1.25, 5.95), arrowprops=dict(arrowstyle='-|>', color='#17456e', lw=1.6))
ax.annotate('', xy=(2.5, 4.9), xytext=(3.75, 5.95), arrowprops=dict(arrowstyle='-|>', color='#17456e', lw=1.6))
ax.add_patch(Rectangle((0.85, 3.4), 3.3, 1.45, fc='white', ec='#2ca02c', lw=1.8))
ax.text(2.5, 4.12, 'recipe (α, β, γ)\ndet JF = bc  ✓ KELLER', ha='center', va='center', fontsize=9.5, color='#1c6b1c')
# right: 2D
ax.add_patch(Rectangle((5.45, 5.6), 4.1, 3.9, fc='#fff0ec', ec='#bb3311', lw=1.6))
ax.text(7.5, 9.05, 'C², weights (1,−1)', ha='center', fontsize=10.5, fontweight='bold', color='#8c2208')
ax.text(7.5, 8.35, 'weight-0 monomials:\n{A−B=0} = {v$^B$}   (one family)', ha='center', fontsize=9.5,
        va='top', bbox=dict(fc='white', ec='#bb3311', boxstyle='round,pad=0.35'))
ax.text(7.5, 6.15, 'v = xy   (only channel)', fontsize=11, ha='center', color='#8c2208')
ax.annotate('', xy=(7.5, 4.9), xytext=(7.5, 5.95), arrowprops=dict(arrowstyle='-|>', color='#8c2208', lw=1.6))
ax.add_patch(Rectangle((5.85, 3.25), 3.3, 1.6, fc='white', ec='#bb3311', lw=1.8))
ax.text(7.5, 4.05, 'shadow (α, β)\ndet = y³ E(v)  ✗ never\na nonzero constant', ha='center', va='center', fontsize=9.5, color='#8c2208')
# cut the channel
ax.add_patch(Rectangle((2.1, 0.5), 5.8, 1.9, fc='#f7f2e7', ec='#8a7350', lw=1.5))
ax.text(5.0, 1.45, 'cut the second channel (b := 0 in γ = 1+av+bt):\ndet JF₃ ≡ 0 exactly  — one channel cannot braid', ha='center', va='center', fontsize=9.8, color='#5e4a2a')

# ---------------- (c) real sign mosaic of the 2-D shadow ----------------
ax = axes[1, 0]; ax.set_facecolor('#fdfcf8')
r0 = -1.43906088216892392
n = 900
g = np.linspace(-3, 3, n)
Xg, Yg = np.meshgrid(g, g)
Vg = Xg * Yg
det_sign = np.sign(Yg**3 * Epoly(Vg))
ax.pcolormesh(Xg, Yg, det_sign, cmap='coolwarm', shading='auto', alpha=0.75)
gg = np.linspace(0.06, 3, 400)
for s in (1, -1):
    yy = gg * s * np.sign(r0)
    xx = r0 / yy
    m = (np.abs(xx) <= 3) & (np.abs(yy) <= 3)
    ax.plot(xx[m], yy[m], 'k-', lw=2.4)
ax.axhline(0, color='k', lw=2.4)
ax.plot(0, 0, 'ko', ms=5)
ax.text(-2.9, 2.55, 'det > 0', color='#b2182b', fontsize=12, fontweight='bold')
ax.text(-2.9, -2.75, 'det < 0', color='#2166ac', fontsize=12, fontweight='bold')
ax.text(0.35, -2.9, f'degeneracy: y = 0 and hyperbola xy = r₀ = {r0:.9f}', fontsize=9.5)
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3)
ax.set_xlabel('x', fontsize=10.5); ax.set_ylabel('y', fontsize=10.5)
ax.set_title('(c) THE SEEDED 2-D SHADOW over ℝ²:\ndet = y³E(xy);  E < 0 for v > r₀;  fiber exactly 4', fontsize=10.5)

# ---------------- (d) Moh's horizon ----------------
ax = axes[1, 1]; ax.set_facecolor('#fdfcf8')
ax.add_patch(Rectangle((0, 0), 100, 10, fc='#e3f3e3', ec='none'))
ax.axvspan(100, 132, color='#f6d9d9')
ax.axvline(100, color='#1c6b1c', lw=2, ls='--')
ax.text(50, 9.15, 'MOH 1983: JC TRUE for degrees \u2264 100  (Crelle 340: 140\u2013212)',
        ha='center', fontsize=9.5, color='#1c6b1c')
ax.text(116, 9.15, 'deg > 100\nOPEN SEA', ha='center', va='center', fontsize=9.5, color='#a11c1c')
# -- tower darts row (y ~ 6..8)
ax.plot([7], [7.9], marker='*', ms=15, color='#d62728', mfc='none')
ax.text(9.5, 7.9, 'fiber-4 tower G:  deg 12, det = \u22126  (3-D)', fontsize=8.2, va='center')
ax.plot([7], [6.5], marker='*', ms=15, color='#d62728')
ax.text(9.5, 6.5, 'fiber-3 tower:  deg 7, det = \u22122  (3-D)', fontsize=8.2, va='center')
wall_x = 63
ax.add_patch(Rectangle((wall_x - 1.3, 5.3), 2.6, 3.5, fc='#c9c9c9', ec='#555', hatch='///'))
ax.annotate('', xy=(wall_x - 1.5, 6.9), xytext=(52, 6.9),
            arrowprops=dict(arrowstyle='-|>', lw=2, color='#d62728'))
ax.text(wall_x, 4.85, 'THE CANAL WALL (this note):\ndet = y\u00b3E(v) \u2260 const', ha='center', va='top', fontsize=8.0, color='#333')
ax.text(wall_x + 3.2, 6.9, 'no tower-shaped\ndescent exists', fontsize=8.6, va='center', ha='left', color='#a11c1c')
# -- Pinchuk row (y ~ 2.8..4.2)
ax.plot([10], [3.9], marker='s', ms=9, color='#7b3294')
ax.plot([25], [3.9], marker='s', ms=9, color='#7b3294')
ax.plot([35], [3.9], marker='s', ms=8, mfc='none', mec='#7b3294')
ax.text(10, 4.35, '10', fontsize=8, ha='center', color='#5b2a72')
ax.text(25, 4.35, '25', fontsize=8, ha='center', color='#5b2a72')
ax.text(35, 3.25, '35', fontsize=8, ha='center', color='#5b2a72')
ax.text(2, 2.55, 'Pinchuk 1994 (real): component degrees (10, 25), total 35 \u2014\nnonvanishing Jacobian, non-injective \u2014 but det NONCONSTANT:\nthe real chamber is dead, the complex Keller chamber stays sealed',
        fontsize=8.2, ha='left', va='top', color='#5b2a72')
# -- avatar row (y ~ 0..1.5)
ax.text(2, 1.35, 'entire avatars live everywhere (Vitushkin-style shears):  F = (e\u02e3, (y + p(e\u00b2\u02e3))e\u207b\u02e3),\ndet \u2261 1, fibers \u2124-lattices, missed set the line {X = 0} \u2014 but not polynomial,\nand not the recipe: the seed is transplanted, not descended',
        fontsize=7.6, ha='left', va='top', color='#225588')
ax.set_xlim(-1, 133); ax.set_ylim(0, 10)
ax.set_xlabel('polynomial degree', fontsize=10.5)
ax.set_yticks([])
ax.set_title('(d) MOH’S HORIZON: the sealed chamber\nand the wall our tower hits', fontsize=10.5)

fig.suptitle('NOTE 18 — THE MISSING BRAID: why the weighted-lift tower has no 2-D shadow, even entire', fontsize=12.5, y=0.995)
fig.tight_layout(rect=[0, 0, 0.99, 0.965])
fig.savefig('/home/user/keller2d_figure.png', dpi=170)
print('saved keller2d_figure.png')
