"""
Note 10 figure (2x2): the generic atlas.
(a) d=3 anchor: wall + strata + real-root regions
(b) generic seed certification table (8/8)
(c) escape slope SPECTRUM: 1/2, 2/3, 3/4 (+1/4 root drift)
(d) the forced swallowtail: wall zoom showing the E6 ramphoid cusp at (1/8, -1/768)
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import sympy as sp

w, s, r = sp.symbols("w s r")
fig, axes = plt.subplots(2, 2, figsize=(12.5, 11.5))

# (a) d=3 anchor
ax = axes[0, 0]
p3 = -w**3 - sp.Rational(3,2)*w**2 + sp.Rational(3,2)*w
D3 = sp.sympify("64*r**3 + 96*r**2*s + 207*r**2 - 114*r*s**2 - 117*r*s + 27*r + 27*s**4 + 35*s**3 - 9*s**2")
D3l = sp.lambdify((s, r), D3, "numpy")
Phi3 = sp.expand(sp.integrate(p3, w))
n = 160
S = np.linspace(-4, 4, n); R = np.linspace(-4, 4, n)
grid = np.zeros((n, n), dtype=int)
for i, rv in enumerate(R):
    for j, sv in enumerate(S):
        hh = sp.Poly(sp.expand((Phi3 - s*w + r).subs({s: float(sv), r: float(rv)})), w)
        roots = np.roots([complex(c) for c in hh.all_coeffs()])
        grid[i, j] = sum(1 for z in roots if abs(z.imag) < 1e-7)
SS, RR = np.meshgrid(S, R)
cmap = ListedColormap(["#2b2b3a", "#dfe8f5", "#b8d8b8", "#90ee90", "#ffd24a"])
ax.pcolormesh(SS, RR, grid, cmap=cmap, shading="auto", vmin=0, vmax=4)
ax.contour(SS, RR, D3l(SS, RR), levels=[0], colors="crimson", linewidths=1.4)
CUSPS3 = [(-2.29904, 1.33702), (0.29904, 0.03798)]
for c in CUSPS3: ax.plot(*c, "*", color="gold", mec="k", ms=17, zorder=17)
ax.plot(-1, -1, "s", color="navy", mec="w", ms=12, zorder=17)
ax.annotate("crunode (-1,-1): contacts t = -2 and t = 1 EXACTLY", xy=(-1,-1), xytext=(-3.9,-2.9),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", color="k"))
legend = [Line2D([0],[0], marker="*", color="w", mfc="gold", mec="k", ms=15, ls="", label="2 cusps (both real)"),
          Line2D([0],[0], marker="s", color="w", mfc="navy", mec="w", ms=11, ls="", label="1 crunode"),
          Line2D([0],[0], color="crimson", lw=1.5, label="wall $D_3=0$ (quartic, 10 terms)")]
ax.legend(handles=legend, loc="upper left", fontsize=8.5, framealpha=0.95)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) the anchor: d=3 chamber (fiber 4) - 2 cusps + 1 node = 3,\nbudget (4-1)(4-2)/2 = 3 balanced; monodromy S4 (|G| = 24, certified)", fontsize=10)

# (b) generic atlas table
ax = axes[0, 1]
ax.axis("off")
rep = json.load(open("generic_atlas.json"))
rows = []
for name, v in rep.items():
    rows.append([name.replace(" (", "\n("), f"{v['wall_deg']}", f"{v['wall_terms']}",
                 f"{v['cusps']} ({v['real_cusps']}R)", f"{v['nodes']} ({v['real_nodes']}R)",
                 "OK" if v["balanced"] else "FAIL", "OK" if v["empty"] else "FAIL"])
cols = ["seed (c3,c4)", "wall\ndeg", "terms", "cusps", "nodes", "budget", "emptiness"]
tab = ax.table(cellText=rows, colLabels=cols, cellLoc="center", loc="center")
tab.auto_set_font_size(False); tab.set_fontsize(9.5); tab.scale(1, 1.9)
for i in range(len(rows)):
    for j in (5, 6):
        tab[(i+1, j)].set_facecolor("#d6f5d6")
for j in range(len(cols)):
    tab[(0, j)].set_facecolor("#dfe8f5")
ax.set_title("(b) GENERIC ATLAS: 8/8 normalized deg-4 seeds pass\nwall deg 5, terms 14, budget 3 cusps + 3 nodes, emptiness - the tower is TYPICAL;\nmonodromy S5 spot-checked on seed (2,-1) (|G| = 120)", fontsize=10)

# (c) escape spectrum
ax = axes[1, 0]
fold = json.load(open("atlas5_escape_fold.json"))
cusp = json.load(open("atlas5_escape_cusp.json"))
swall = json.load(open("swallowtail_escape.json"))
drift = json.load(open("swallowtail_drift.json"))
ax.loglog(fold["delta"], fold["gamma"], "o-", color="navy", ms=4, lw=1.2,
          label=f"fold m=2: slope {fold['slope']:.3f} ~ 1/2")
ax.loglog(cusp["delta"], cusp["gamma"], "s-", color="crimson", ms=4, lw=1.2,
          label=f"cusp m=3: slope {cusp['slope']:.3f} ~ 2/3")
ax.loglog(swall["delt"], swall["gam"], "D-", color="darkgreen", ms=5, lw=1.2,
          label=f"SWALLOWTAIL m=4 (forced seed): slope {swall['slope']:.3f} ~ 3/4")
ax.loglog(drift["delt"], drift["drift"], "^:", color="purple", ms=5, lw=1.2,
          label=f"root drift |w-t0|: slope {drift['slope']:.3f} ~ 1/4 = 1/m")
dref = np.array(fold["delta"])
for k_, c_ in [(0.5, "navy"), (2/3, "crimson"), (0.75, "darkgreen")]:
    const = {0.5: 0.5, 2/3: 1.2, 0.75: 0.87}[k_]
    ax.loglog(dref, const*dref**k_, "--", color=c_, alpha=0.35, lw=1)
ax.set_xlabel("distance delta to wall"); ax.set_ylabel("|gamma|  or  |w - t0|")
ax.legend(fontsize=9, loc="lower right")
ax.set_title("(c) THE ESCAPE LAW: slope = (m-1)/m for a multiplicity-m fiber contact\n"
             "1/2 and 2/3 in every tower chamber; 3/4 forced - root drift confirms 1/m", fontsize=10)

# (d) swallowtail wall zoom
ax = axes[1, 1]
DSW = sp.sympify(open("swallowtail_wall.txt").read())
DSWl = sp.lambdify((s, r), DSW, "numpy")
xs = np.linspace(0.04, 0.24, 700); ys = np.linspace(-0.03, 0.03, 500)
SS2, RR2 = np.meshgrid(xs, ys)
DD2 = DSWl(SS2, RR2)
ax.contour(SS2, RR2, DD2, levels=[0], colors="crimson", linewidths=1.3)
ax.plot(0.125, -1/768, "*", color="purple", mec="k", ms=22, zorder=17)
ax.plot(0.156416, 0.00860556, "*", color="gold", mec="k", ms=16, zorder=17)
ax.annotate("A2 cusp (ordinary)", xy=(0.156416, 0.00860556), xytext=(0.175, 0.022),
            fontsize=9, arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("E6 RAMPHOID CUSP (swallowtail):\nh = (2w-1)^4 (8w^2-104w-1)/768,\n"
            "fiber mult 4, wall semigroup <3,4>, delta = 3", xy=(0.125, -1/768),
            xytext=(0.045, -0.0275), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="k"))
ax.set_xlabel("s"); ax.set_ylabel("r")
ax.set_title("(d) the forced swallowtail in its wall (zoom): the violet star is a\ngenuine non-empty (4,1,1) stratum - budget 3 + 2 + 5 = 10 still balances", fontsize=10)

fig.suptitle("The generic atlas: the tower is typical - and a forced swallowtail nails the (m-1)/m escape law", fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig("generic_atlas.png", dpi=150)
print("saved generic_atlas.png")
