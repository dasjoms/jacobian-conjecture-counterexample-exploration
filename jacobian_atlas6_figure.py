"""
Note 8 figure (2x2):
(a) real region atlas {1,3} real roots; wall D6=0; real cusp (whisker), 2 acnodes
(b) cusp loop braid (3-cycle bold via actual permutation)
(c) escape rates: fold 0.4988, cusp 0.6635 (universality)
(d) real census: {1: 82.54%, 3: 17.46%} - odd only, nothing missed
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import sympy as sp

s_, r_ = sp.symbols("s r")
D6 = sp.sympify(open("atlas6_wall.txt").read())
disc = sp.lambdify((s_, r_), D6, "numpy")

fig, axes = plt.subplots(2, 2, figsize=(12.5, 11.5))

# (a) region atlas
ax = axes[0, 0]
gz = np.load("atlas6_grid.npz"); S, R, grid = gz["S"], gz["R"], gz["grid"]
SS, RR = np.meshgrid(S, R)
cmap = ListedColormap(["#2b2b3a", "#dfe8f5", "#ffd24a", "#b8d8b8", "#f5c6f5"])
ax.pcolormesh(SS, RR, grid, cmap=cmap, shading="auto", vmin=0, vmax=4)
DD = disc(SS, RR)
ax.contour(SS, RR, DD, levels=[0], colors="crimson", linewidths=1.4)
CUSP = (0.3043412601764495, 0.03338254261374218)
ACN = [(-1.76416996, 2.29113957), (-0.10068230, -0.86419633)]
ax.plot(*CUSP, "*", color="gold", mec="k", ms=19, zorder=17)
for an in ACN:
    ax.plot(*an, "D", color="darkviolet", mec="w", ms=10, zorder=17)
ax.annotate("whisker root:\ncusp target missed over R", xy=CUSP, xytext=(0.9, 1.3),
            fontsize=9, arrowprops=dict(arrowstyle="->", color="k"))
legend = [Line2D([0],[0], marker="*", color="w", mfc="gold", mec="k", ms=16, ls="", label="real cusp (missed over R: 0/4 residual roots real)"),
          Line2D([0],[0], marker="D", color="w", mfc="darkviolet", mec="w", ms=9, ls="", label="acnodes (2), escaped via complex contacts"),
          Line2D([0],[0], color="crimson", lw=1.5, label="wall $D_6=0$ (irreducible septic)")]
ax.legend(handles=legend, loc="upper left", fontsize=8.5, framealpha=0.95)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) real fiber atlas: regions by #real preimages (1=blue, 3=gold)\n"
             "NO 0-region: odd fiber degree => no missed cone; strata: 1 real cusp, 2 acnodes", fontsize=10)

# (b) cusp braid
ax = axes[0, 1]
steps = 2400
ts = np.linspace(0, 2*np.pi, steps)
path = [(CUSP[0] + 0.006*np.exp(1j*t)/np.sqrt(2), CUSP[1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)) for t in ts]
def hc(sv, rv): return [-1/7, 1/6, 0, 0, -20/21, 13/14, -sv, rv]
def greedy(prev, nxt):
    used = [False]*7; out = [None]*7
    for a in range(7):
        bd, bj = 1e30, 0
        for b in range(7):
            if not used[b] and abs(prev[a]-nxt[b]) < bd:
                bd, bj = abs(prev[a]-nxt[b]), b
        used[bj] = True; out[a] = nxt[bj]
    return out
traces = [list(np.roots(hc(*path[0])))]
for i in range(1, steps):
    traces.append(greedy(traces[-1], list(np.roots(hc(*path[i])))))
init = traces[0]
used = [False]*7; perm = [None]*7
for a in range(7):
    bd, bj = 1e30, 0
    for b in range(7):
        if used[b]: continue
        d = abs(traces[-1][a]-init[b])
        if d < bd: bd, bj = d, b
    used[bj] = True; perm[a] = bj
moving = set()
for i in range(7):
    if perm[i] != i:
        j = i
        while j not in moving:
            moving.add(j); j = perm[j]
print("cusp loop perm:", perm, " moving:", sorted(moving))
traces.append(traces[0])
hot_cols = ["crimson", "navy", "darkgreen"]
for k in range(7):
    zs = [tr[k] for tr in traces]
    hot = k in moving
    ax.plot([z.real for z in zs], [z.imag for z in zs], "-",
            color=hot_cols[sorted(moving).index(k)] if hot else "#999999",
            lw=2.6 if hot else 0.8, alpha=0.95 if hot else 0.5)
for k in range(7):
    ax.plot(init[k].real, init[k].imag, "o", color="k", ms=5)
ax.set_xlabel("Re w"); ax.set_ylabel("Im w")
ax.set_title("(b) cusp loop braid: 3 sheets cycle (bold), 4 stay put\nmonodromy = (5 7 6), S7-certified", fontsize=10)

# (c) escape rates
ax = axes[1, 0]
fold = json.load(open("atlas6_escape_fold.json")); cusp = json.load(open("atlas6_escape_cusp.json"))
ax.loglog(fold["delta"], fold["gamma"], "o-", color="navy", ms=4, lw=1,
          label=f"fold: slope {fold['slope']:.3f} (~1/2), |gamma| ~ {fold['const']:.3f} delta^0.50")
ax.loglog(cusp["delta"], cusp["gamma"], "s-", color="crimson", ms=4, lw=1,
          label=f"cusp: slope {cusp['slope']:.3f} (~2/3), |gamma| ~ {cusp['const']:.3f} delta^0.66")
dref = np.array(fold["delta"])
ax.loglog(dref, 0.52*dref**0.5, "--", color="navy", alpha=0.4, lw=1)
ax.loglog(dref, 1.23*dref**(2/3), "--", color="crimson", alpha=0.4, lw=1)
ax.set_xlabel("distance delta to wall"); ax.set_ylabel("|gamma| on escaping sheet (|x| = |C|/|gamma|)")
ax.legend(fontsize=9, loc="lower right")
ax.set_title("(c) escape universality across chambers: fold delta^{-1/2}, cusp delta^{-2/3}\n(F4: 1/2; F5: 0.4986/0.6635; F6: 0.4988/0.6635)", fontsize=10)

# (d) census
ax = axes[1, 1]
cen = json.load(open("atlas6_realcensus.json"))
keys = ["0","1","2","3","4","5","6","7"]
vals = [100*cen.get(k,0)/200000 for k in keys]
cols = ["#2b2b3a", "#4a7ebb", "#aaaaaa", "#ffd24a", "#aaaaaa", "#b8d8b8", "#aaaaaa", "#f5c6f5"]
bars = ax.bar(keys, vals, color=cols, ec="k")
for k, v, b in zip(keys, vals, bars):
    ax.text(b.get_x()+b.get_width()/2, v+0.7, f"{v:.2f}%" if v else "0", ha="center", fontsize=9)
ax.set_xlabel("# of REAL preimages of a random real target (A,B,1)")
ax.set_ylabel("percent of 200,000 targets")
ax.set_ylim(0, 96)
ax.set_title("(d) R-side reality in an ODD chamber: every target hit (0.0000% missed);\n"
             "only odd counts occur; 5 or 7 preimages never drawn (rare/absent regions)", fontsize=10)

fig.suptitle("The septic chamber: wall of $F_6$ - 5 cusps + 10 nodes = 15, braided by $S_7$", fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig("septic_atlas.png", dpi=150)
print("saved septic_atlas.png")
