"""
Note 7 figure (2x2):
(a) real (s,r) atlas: regions colored by #real roots of h; wall overlay; strata marked
(b) cusp loop braid: root traces around the cusp (0.29716, 0.03303) -> 3-cycle
(c) escape rates: fold (slope 1/2) vs cusp (slope 2/3), log-log
(d) real census: histogram {0: 8.69%, 2: 84.91%, 4: 6.39%}, 6 never observed
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

D5_terms = [(3037500000,5,0), (-3037500000,4,1), (11907000000,4,0), (-50625000,3,2),
            (-19265850000,3,1), (15615855000,3,0), (9797625000,2,3), (8794980000,2,2),
            (-23436459000,2,1), (8496467856,2,0), (-6968981250,1,4), (2780676000,1,3),
            (8071107300,1,2), (-6147920736,1,1), (892142910,1,0), (1220703125,0,6),
            (-46159500,0,5), (-2098753050,0,4), (1564980260,0,3), (-247817475,0,2)]
def disc_val(s, r):
    return sum(c * r**rr * s**ss for c, rr, ss in D5_terms)

fig, axes = plt.subplots(2, 2, figsize=(12.5, 11.5))

# ---- (a) real atlas ----
ax = axes[0, 0]
gz = np.load("atlas5_grid.npz")
S, R, grid = gz["S"], gz["R"], gz["grid"]
SS, RR = np.meshgrid(S, R)
cmap = ListedColormap(["#2b2b3a", "#dfe8f5", "#b8d8b8", "#ffe9a8"])
ax.pcolormesh(SS, RR, grid, cmap=cmap, shading="auto")
DD = disc_val(SS, RR)
ax.contour(SS, RR, DD, levels=[0], colors="crimson", linewidths=1.4)
strata = {
    "cusp (real)": [(0.2971583157757676, 0.03302630284803969), (-2.6520401043053803, 1.18422351654555)],
    "crunode": [(-0.8818834118879952, -0.8833717438605823)],
    "acnode": [(0.4920279747898448, -0.4640637734449852)],
}
for (sx, rx) in strata["cusp (real)"]:
    ax.plot(sx, rx, "*", color="gold", mec="k", ms=17, zorder=17)
ax.plot(*strata["crunode"][0], "s", color="navy", mec="w", ms=11, zorder=17)
ax.plot(*strata["acnode"][0], "D", color="darkviolet", mec="w", ms=10, zorder=17)
from matplotlib.lines import Line2D
legend = [Line2D([0],[0], marker="*", color="w", mfc="gold", mec="k", ms=15, ls="", label="cusps (2 real)"),
          Line2D([0],[0], marker="s", color="w", mfc="navy", mec="w", ms=10, ls="", label="crunode (real crossings)"),
          Line2D([0],[0], marker="D", color="w", mfc="darkviolet", mec="w", ms=9, ls="", label="acnode (isolated real point)"),
          Line2D([0],[0], color="crimson", lw=1.5, label="wall $D_5=0$")]
ax.legend(handles=legend, loc="upper left", fontsize=9, framealpha=0.95)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) real fiber atlas: regions by #real preimages (0=dark, 2=blue, 4=green)\n"
             "6-real region absent in window; wall has 2 cusps, 1 crunode, 1 acnode", fontsize=10)
for (sx, rx), num in [((-2.6520401043, 1.1842235165), None)]:
    pass
# annotate region counts at sample points
for (sx, rx, lab) in [(2.6, 2.6, "2"), (0.5, 2.6, "2"), (-3.0,-2.5,"0"), (1.5,-0.35,"0?"), (3.2,-0.2,"4"), (-2.2,-0.5,"4")]:
    pass

# ---- (b) cusp braid ----
ax = axes[0, 1]
cu = (0.2971583157757676, 0.03302630284803969)
steps = 2400
ts = np.linspace(0, 2*np.pi, steps)
path = [(cu[0] + 0.008*np.exp(1j*t)/np.sqrt(2), cu[1] + 1j*0.008*np.exp(1j*t)/np.sqrt(2)) for t in ts]
def greedy(prev, nxt):
    used = [False]*6; out = [None]*6
    for a in range(6):
        bd, bj = 1e30, 0
        for b in range(6):
            if not used[b] and abs(prev[a]-nxt[b]) < bd:
                bd, bj = abs(prev[a]-nxt[b]), b
        used[bj] = True; out[a] = nxt[bj]
    return out
traces = [list(np.roots([-1/6, 1/5, 0, -14/15, 9/10, -path[0][0], path[0][1]]))]
for i in range(1, steps):
    traces.append(greedy(traces[-1], list(np.roots([-1/6, 1/5, 0, -14/15, 9/10, -path[i][0], path[i][1]]))))
# identify the 3-cycle indices from the end-to-start matching permutation
init = traces[0]
used = [False]*6; perm = [None]*6
for a in range(6):
    bd, bj = 1e30, 0
    for b in range(6):
        if used[b]: continue
        d = abs(traces[-1][a]-init[b])
        if d < bd: bd, bj = d, b
    used[bj] = True; perm[a] = bj
moving = set()
for i in range(6):
    if perm[i] != i:
        j = i
        while j not in moving:
            moving.add(j); j = perm[j]
print("perm:", perm, " moving sheets:", sorted(moving))
traces.append(traces[0])  # close
cols = {3: "crimson", 4: "navy", 5: "darkgreen", 0: "crimson", 1: "navy", 2: "darkgreen"}
for k in range(6):
    zs = [tr[k] for tr in traces]
    hot = k in moving
    ax.plot([z.real for z in zs], [z.imag for z in zs], "-",
            color=cols[sorted(moving).index(k)] if hot else "#999999",
            lw=2.6 if hot else 0.8, alpha=0.95 if hot else 0.5)
z0 = traces[0]
for k in range(6):
    ax.plot(z0[k].real, z0[k].imag, "o", color="k", ms=5)
ax.set_xlabel("Re w"); ax.set_ylabel("Im w")
ax.set_title("(b) cusp loop braid: 3 sheets cycle (bold), 3 stay put\nmonodromy = (4 6 5), a 3-cycle (like note 6)", fontsize=10)

# ---- (c) escape rates ----
ax = axes[1, 0]
fold = json.load(open("atlas5_escape_fold.json")); cusp = json.load(open("atlas5_escape_cusp.json"))
ax.loglog(fold["delta"], fold["gamma"], "o-", color="navy", ms=4, lw=1,
          label=f"fold: slope {fold['slope']:.3f} (~1/2), |gamma|~0.50 delta^0.50")
ax.loglog(cusp["delta"], cusp["gamma"], "s-", color="crimson", ms=4, lw=1,
          label=f"cusp: slope {cusp['slope']:.3f} (~2/3), |gamma|~1.19 delta^0.66")
dref = np.array(fold["delta"])
ax.loglog(dref, 0.5*dref**0.5, "--", color="navy", alpha=0.4, lw=1)
ax.loglog(dref, 1.2*dref**(2/3), "--", color="crimson", alpha=0.4, lw=1)
ax.set_xlabel("distance delta to wall"); ax.set_ylabel("|gamma| on escaping sheet (|x| = |C|/|gamma|)")
ax.legend(fontsize=9, loc="lower right")
ax.set_title("(c) escape rates: |x| blows up like delta^(-1/2) at folds,\nlike delta^(-2/3) at cusps [new exponent for the sextic]", fontsize=10)

# ---- (d) real census ----
ax = axes[1, 1]
cen = json.load(open("atlas5_realcensus.json"))
keys = ["0","1","2","3","4","5","6"]
vals = [100*cen.get(k,0)/200000 for k in keys]
bars = ax.bar(keys, vals, color=["#2b2b3a","#aaaaaa","#4a7ebb","#aaaaaa","#7eb77e","#aaaaaa","#ffd24a"], ec="k")
for k, v, b in zip(keys, vals, bars):
    ax.text(b.get_x()+b.get_width()/2, v+0.6, f"{v:.2f}%", ha="center", fontsize=10)
ax.set_xlabel("# of REAL preimages of a random real target (A,B,1)")
ax.set_ylabel("percent of 200,000 targets")
ax.set_ylim(0, 97)
ax.set_title("(d) R-side reality: ~8.7% of real targets have NO real preimage (missed!);\nodd counts never occur; 6 preimages never observed in 200k draws", fontsize=10)

fig.suptitle("The escape atlas of $F_5$: a maximally singular sextic wall, braided by $S_6$", fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig("sextic_atlas.png", dpi=150)
print("saved sextic_atlas.png")
