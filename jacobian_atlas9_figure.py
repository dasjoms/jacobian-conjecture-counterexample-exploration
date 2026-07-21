"""
Note 11 figure: the decic chamber atlas.
(a) real atlas F9 in the (s,r) plane: region grid (real-fiber count), wall curve,
    2 REAL CUSPS (both hit), crunode = missed cone corner, 3 acnodes.
(b) wall monomial support vs the weighted cone (n-1)j + n i <= n(n-1), n=10:
    56 cone sites, 53 occupied, THREE HOLES marked: (0,0), (0,1), (8,0).
(c) real census across the five even-n chambers surveyed (n=4,6,8,10) + odd (n=9):
    missed-region mass ~8.7% in every even chamber.
(d) prediction-ledger strip (text).
"""
import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

fig = plt.figure(figsize=(15.5, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.34, wspace=0.22)

# ---------- (a) real atlas ----------
ax = fig.add_subplot(gs[0, 0])
g = np.load("atlas9_grid.npz")
S, R, G = g["S"], g["R"], g["grid"]
ax.imshow(G, extent=[S.min(), S.max(), R.min(), R.max()], origin="lower",
          aspect="auto", cmap="viridis", alpha=0.85, vmin=0, vmax=4,
          interpolation="nearest")
# wall curve: contour D9 = 0 on the same grid via marching squares on sign of D9
import sympy as sp
sr, rr = sp.symbols("s r")
D9 = sp.sympify(open("atlas9_wall.txt").read())
fD = sp.lambdify((sr, rr), D9, "numpy")
SS, RR = np.meshgrid(np.linspace(S.min(), S.max(), 700), np.linspace(R.min(), R.max(), 700))
with np.errstate(all="ignore"):
    Z = np.sign(fD(SS, RR).astype(np.float64))
ax.contour(SS, RR, Z, levels=[0.0], colors=["#ff2233"], linewidths=[2.2])
# singularities
_cx = json.load(open("atlas9_cuspfix.json"))
cusp1 = (float(_cx["d9_cusp1"]["s0"]), float(_cx["d9_cusp1"]["r0"]))   # t = +0.3299 (hit cusp)
cusp2 = (float(_cx["d9_cusp2"]["s0"]), float(_cx["d9_cusp2"]["r0"]))   # t = -0.8914 (hit cusp)
crun = (-0.928716289105, -0.929242511862)
acn = [(-1.46227202662, 2.16224722344), (1.12101982377, 0.449304050481), (-0.764724180785, -1.09356697414)]
for (x, y), lab in [(cusp1, "cusp (hit)"), (cusp2, "cusp (hit)")]:
    ax.plot([x], [y], marker="*", ms=17, mfc="gold", mec="k", zorder=6)
for x, y in acn:
    ax.plot([x], [y], marker="o", ms=8, mfc="none", mec="#00e5ff", mew=2.2, zorder=6)
ax.plot([crun[0]], [crun[1]], marker="X", ms=13, mfc="#ff9100", mec="k", zorder=7)
ax.annotate("crunode = missed cone corner\n(residual: 0 real roots)", crun,
            textcoords="offset points", xytext=(12, -34), fontsize=9, color="#ff9100",
            arrowprops=dict(arrowstyle="->", color="#ff9100"))
ax.plot([], [], "*", ms=13, mfc="gold", mec="k", label="real cusps (2, both HIT)")
ax.plot([], [], "o", ms=8, mfc="none", mec="#00e5ff", mew=2.2, label="acnodes (3)")
ax.plot([], [], "X", ms=11, mfc="#ff9100", mec="k", label="crunode (1)")
ax.plot([], [], "-", color="#ff2233", lw=2.2, label="wall D9 (deg 10, 53 terms)")
ax.set_title("(a) F9 real atlas: 0/2/4-real-fiber regions, wall, and its 6 real singular points",
             fontsize=11)
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
ax.set_xlabel("s"); ax.set_ylabel("r")

# ---------- (b) support map ----------
ax = fig.add_subplot(gs[0, 1])
n = 10
P = sp.Poly(D9, sr, rr)
have = {(int(m[0][0]), int(m[0][1])) for m in P.terms()}   # (j=s-exp, i=r-exp) -> plot x=j,y=i
# LEDGER: an earlier draft transposed this set; the r^8 hole rendered at the wrong corner.
for i in range(n + 1):
    for j in range(n + 1):
        if (n - 1) * j + n * i <= n * (n - 1):
            occ = (j, i) in have
            hole = not occ
            ax.plot([j], [i], "s", ms=9,
                    mfc=("#1a7f37" if occ else "white"),
                    mec=("#1a7f37" if occ else "#d32f2f"), mew=1.8, zorder=4 if occ else 5)
# frontier line
jj = np.linspace(0, n, 200)
ax.plot(jj, (n * (n - 1) - (n - 1) * jj) / n, "k--", lw=1.4, alpha=0.65,
        label=r"cone edge $(n\!-\!1)j + n\,i = n(n\!-\!1)$")
for (j, i), lab, off in [((0, 0), "(const)", (12, -13)), ((1, 0), "($s^1$)", (8, 10)), ((0, 8), "($r^{8}$)", (12, -10))]:
    ax.annotate(lab, (j, i), textcoords="offset points", xytext=off,
                fontsize=10, color="#d32f2f", fontweight="bold")
ax.plot([], [], "s", ms=9, mfc="#1a7f37", mec="#1a7f37", label="present (53 monomials)")
ax.plot([], [], "s", ms=9, mfc="white", mec="#d32f2f", mew=1.8, label="HOLE (3)")
ax.set_xlabel("j  (power of s)"); ax.set_ylabel("i  (power of r)")
ax.set_title("(b) wall D9 support vs weight cone: 56 cone sites − 3 holes = 53 terms", fontsize=11)
ax.legend(loc="upper right", fontsize=8.5); ax.set_xlim(-0.6, 10.8); ax.set_ylim(-0.6, 9.8)

# ---------- (c) census panorama ----------
ax = fig.add_subplot(gs[1, 0])
# census values from notes 9-11 (200k samples each, normal(0,1.5^2), C=1)
def pct(fn, k):
    d = json.load(open(fn)); return 100.0 * d[str(k)] / sum(d.values())
cens = [
    ("n=6 (d5)",  None, {0: 8.69, 2: 84.91, 4: 6.39}),
    ("n=7 (d6)",  None, {1: 82.54, 3: 17.46}),
    ("n=8 (d7)",  None, {0: 8.74, 2: 84.27, 4: 7.00}),
    ("n=9 (d8)",  None, {1: 82.60, 3: 17.40}),
    ("n=10 (d9)", "atlas9_realcensus.json", None),
]
xs, missed, other1, other2, labs = [], [], [], [], []
for name, fn, manual in cens:
    if fn:
        m = pct(fn, 0); vals = {0: m, 2: pct(fn, 2), 4: pct(fn, 4)}
    else:
        vals = manual
    ks = sorted(vals)
    xs.append(name); missed.append(vals[ks[0]])
    other1.append(vals[ks[1]] if len(ks) > 1 else 0)
    other2.append(vals[ks[2]] if len(ks) > 2 else 0)
xp = np.arange(len(xs))
ax.bar(xp, missed, 0.55, color="#d32f2f", label="0-real (even n) / 1-real (odd n): missed+whisker")
ax.bar(xp, other1, 0.55, bottom=missed, color="#388e3c", label="middle fiber count")
ax.bar(xp, other2, 0.55, bottom=np.array(missed)+np.array(other1), color="#90a4ae", label="top fiber count")
for i, m in enumerate(missed):
    ax.annotate(f"{m:.2f}%", (xp[i], m/2), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
ax.axhline(8.7, color="#d32f2f", ls=":", lw=1)
ax.set_xticks(xp); ax.set_xticklabels(xs, fontsize=9)
ax.set_ylabel("% of sampled targets")
ax.set_title("(c) real census across chambers: the even-n missed cone holds ~8.7% everywhere", fontsize=11)
ax.legend(fontsize=8.5, loc="upper right")

# ---------- (d) prediction ledger strip ----------
ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
rows = [
    ("PREDICTIONS LOCKED BEFORE COMPUTATION (this round)", "k", True),
    ("beta_tower(10)=0  =>  third hole r^8 => terms(10)=53", "g", False),
    ("K(d=9) = den^2 = 225  (integer normalization)", "g", False),
    ("Sturm exact real cusp count d=9: 2", "g", False),
    ("cusp fiber: 3 escaping (1 real) + 7 bounded (1 real)", "g", False),
    ("cusp escape exponent 2/3, all directions, d=7/8/9 (A3 universal)", "g", False),
    ("frozen-point artifact: gamma-floor ~ delta0^(2/3) = 2.6e-4 (obs 2.4e-4)", "g", False),
    ("det JF = 1 exactly at 5/5 rational points: d = 6, 7, 9", "g", False),
    ("monodromy = S10 (Jordan cert: primitive + transposition)", "g", False),
    ("census d=9: missed cone 8.726% (predict 8.6-8.9)", "g", False),
    ("", "k", False),
    ("HONESTY LEDGER (new entries)", "k", True),
    ("hardcoded cusp from 5-dp printout -> bogus slope 0.2746, fiber 10/7", "r", False),
    ("monic-vs-integer normalization printed fake K=False for d=8,9", "r", False),
    ("off-by-one audit rows (h'=p-s0, h''=p', h'''=p'') caught on rerun", "r", False),
]
y = 0.98
for txt, c, bold in rows:
    ax.text(0.02, y, txt, fontsize=10 if bold else 9.2, fontweight="bold" if bold else "normal",
            color={"g": "#1b5e20", "r": "#b71c1c", "k": "black"}[c], va="top", transform=ax.transAxes)
    y -= 0.072
ax.set_title("(d) scoreboard strip (green = prediction redeemed, red = bug caught + fixed)", fontsize=11)

fig.suptitle("JACOBIAN LAB — Note 11: the DECIC chamber (fiber 10) — wall 53 terms, 28 nodes, S10, and the three-hole law",
             fontsize=13.5, fontweight="bold")
fig.savefig("jacobian_atlas9_atlas.png", dpi=150, bbox_inches="tight")
print("saved jacobian_atlas9_atlas.png")
