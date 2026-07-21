"""
Note 9 figure (2x2): octic chamber atlas.
(a) region atlas {0,2,4}; wall; strata: 2 real cusps, 1 crunode (MISSED over R), 2 acnodes
(b) cusp braid, (c) escape slopes, (d) census {0,2,4}
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import sympy as sp

s_, r_ = sp.symbols("s r")
D7 = sp.sympify(open("atlas7_wall.txt").read())
disc = sp.lambdify((s_, r_), D7, "numpy")

fig, axes = plt.subplots(2, 2, figsize=(12.5, 11.5))

ax = axes[0, 0]
gz = np.load("atlas7_grid.npz"); S, R, grid = gz["S"], gz["R"], gz["grid"]
SS, RR = np.meshgrid(S, R)
cmap = ListedColormap(["#2b2b3a", "#f5c6f5", "#dfe8f5", "#f5c25c", "#b8d8b8", "#f5c25c", "#f0b0ff", "#f5c25c", "#90ee90"])
ax.pcolormesh(SS, RR, grid, cmap=cmap, shading="auto", vmin=0, vmax=8)
DD = disc(SS, RR)
ax.contour(SS, RR, DD, levels=[0], colors="crimson", linewidths=1.3)
CUSPS = [(0.3104763509023746, 0.034001356659445606), (-3.0373941, 1.3862193)]
CRU = (-0.90945711, -0.91031417)
ACN = [(0.13989048, 1.54077491), (-0.44418202, -1.01180283)]
for c in CUSPS: ax.plot(*c, "*", color="gold", mec="k", ms=18, zorder=17)
ax.plot(*CRU, "s", color="navy", mec="w", ms=12, zorder=17)
for an in ACN: ax.plot(*an, "D", color="darkviolet", mec="w", ms=10, zorder=17)
ax.annotate("crunode: MISSED over R\n(both real contacts escape,\nresidual quartic 0 real)", xy=CRU, xytext=(-3.6, -2.9),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", color="k"))
legend = [Line2D([0],[0], marker="*", color="w", mfc="gold", mec="k", ms=15, ls="", label="real cusps (2; hit over R - odd residual)"),
          Line2D([0],[0], marker="s", color="w", mfc="navy", mec="w", ms=11, ls="", label="crunode (MISSED over R!)"),
          Line2D([0],[0], marker="D", color="w", mfc="darkviolet", mec="w", ms=9, ls="", label="acnodes (2)"),
          Line2D([0],[0], color="crimson", lw=1.5, label="wall $D_7=0$ (irreducible octic)")]
ax.legend(handles=legend, loc="upper right", fontsize=8, framealpha=0.95)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) real fiber atlas (0=dark missed-cone, 2=blue, 4=green)\n"
             "the cone returns in the even octic chamber; strata: 2 cusps, 1 crunode, 2 acnodes", fontsize=10)

ax = axes[0, 1]
steps = 2400
ts = np.linspace(0, 2*np.pi, steps)
path = [(CUSPS[0][0] + 0.006*np.exp(1j*t)/np.sqrt(2), CUSPS[0][1] + 1j*0.006*np.exp(1j*t)/np.sqrt(2)) for t in ts]
def greedy(prev, nxt):
    used = [False]*8; out = [None]*8
    for a in range(8):
        bd, bj = 1e30, -1
        for b in range(8):
            if not used[b] and abs(prev[a]-nxt[b]) < bd:
                bd, bj = abs(prev[a]-nxt[b]), b
        used[bj] = True; out[a] = nxt[bj]
    return out
def hc(sv, rv): return [-1/8, 1/7, 0, 0, 0, -27/28, 53/56, -sv, rv]
traces = [list(np.roots(hc(*path[0])))]
for i in range(1, steps):
    traces.append(greedy(traces[-1], list(np.roots(hc(*path[i])))))
init = traces[0]
used = [False]*8; perm = [None]*8
for a in range(8):
    bd, bj = 1e30, -1
    for b in range(8):
        if used[b]: continue
        d = abs(traces[-1][a]-init[b])
        if d < bd: bd, bj = d, b
    used[bj] = True; perm[a] = bj
moving = set()
for i in range(8):
    if perm[i] != i:
        j = i
        while j not in moving:
            moving.add(j); j = perm[j]
traces.append(traces[0])
hot_cols = ["crimson", "navy", "darkgreen"]
for k in range(8):
    zs = [tr[k] for tr in traces]
    hot = k in moving
    ax.plot([z.real for z in zs], [z.imag for z in zs], "-",
            color=hot_cols[sorted(moving).index(k)] if hot else "#999999",
            lw=2.6 if hot else 0.8, alpha=0.95 if hot else 0.5)
for k in range(8):
    ax.plot(init[k].real, init[k].imag, "o", color="k", ms=5)
ax.set_xlabel("Re w"); ax.set_ylabel("Im w")
ax.set_title("(b) cusp loop braid: 3-cycle (bold) among 8 sheets\nmonodromy = (6 8 7), S8-certified", fontsize=10)

ax = axes[1, 0]
fold = json.load(open("atlas7_escape_fold.json")); cusp = json.load(open("atlas7_escape_cusp.json"))
ax.loglog(fold["delta"], fold["gamma"], "o-", color="navy", ms=4, lw=1,
          label=f"fold: slope {fold['slope']:.3f}, |gamma| ~ {fold['const']:.3f} delta^0.50")
ax.loglog(cusp["delta"], cusp["gamma"], "s-", color="crimson", ms=4, lw=1,
          label=f"cusp: slope {cusp['slope']:.3f}, |gamma| ~ {cusp['const']:.3f} delta^0.66")
ax.legend(fontsize=9, loc="lower right")
ax.set_xlabel("distance delta to wall"); ax.set_ylabel("|gamma| (|x| = |C|/|gamma|)")
ax.set_title("(c) escape universality, 4th chamber: fold 0.499 ~ 1/2, cusp 0.664 ~ 2/3", fontsize=10)

ax = axes[1, 1]
cen = json.load(open("atlas7_realcensus.json"))
keys = ["0","1","2","3","4","5","6","7","8"]
vals = [100*cen.get(k,0)/200000 for k in keys]
cols = ["#2b2b3a", "#aaaaaa", "#4a7ebb", "#aaaaaa", "#b8d8b8", "#aaaaaa", "#ffd24a", "#aaaaaa", "#f0b0ff"]
bars = ax.bar(keys, vals, color=cols, ec="k")
for k, v, b in zip(keys, vals, bars):
    ax.text(b.get_x()+b.get_width()/2, v+0.7, f"{v:.2f}%" if v else "0", ha="center", fontsize=8.5)
ax.set_xlabel("# of REAL preimages of a random real target (A,B,1)")
ax.set_ylabel("percent of 200,000 targets"); ax.set_ylim(0, 97)
ax.set_title("(d) R-side in the octic chamber: 8.74% missed (cone), even counts only;\n"
             "mirror-image of the sextic census (8.69%)", fontsize=10)

fig.suptitle("The octic chamber: wall of $F_7$ - 6 cusps + 15 nodes = 21, braided by $S_8$; a node-whisker appears", fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig("octic_atlas.png", dpi=150)
print("saved octic_atlas.png")
