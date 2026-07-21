"""
Note 10, stage G: the NONIC ATLAS figure. 4 panels:
(a) real (s,r) atlas: 9-root sea / 3-vs-1 shading, wall, whisker star, acnodes, fold
(b) fiber-count census across the four previous chambers (d=5..8), odd/even rhyme
(c) monodromy certificate: four loop perms + S9 badge + 36/36 transpositions
(d) term law: cone C(n) = n(n+1)/2 + 1, observed terms, holes rule; s^n fingerprints
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# (a) real atlas
dat = np.load("atlas8_grid.npz"); S, R_, grid = dat["S"], dat["R"], dat["grid"]
def p8w(w): return -w**8 + w**7 - 35*w**2/12 + 23*w/12
def Phi9w(w): return -w**9/9 + w**8/8 - 35*w**3/36 + 23*w**2/24
fig, axs = plt.subplots(2, 2, figsize=(15.5, 11.5))
ax = axs[0][0]
ax.pcolormesh(S, R_, np.where(grid==1, 0, np.where(grid==3, 1, grid)), cmap="BuGn",
              shading="auto", alpha=0.55, rasterized=True)
t = np.linspace(-1.7, 1.75, 4000)
sw = p8w(t); rw = t*sw - Phi9w(t)
m = (np.abs(sw) < 4) & (np.abs(rw) < 4)
# plot contiguous segments
start = None
for i in range(len(t)):
    if m[i] and start is None: start = i
    if (not m[i] or i == len(t)-1) and start is not None:
        e = i
        ax.plot(sw[start:e], rw[start:e], "crimson", lw=1.4)
        start = None
CUSP = (0.3151611796143863, 0.034567015896823955)
ax.plot(*CUSP, "*", ms=22, mec="k", mfc="gold", mew=1.1, zorder=5, label="whisker cusp (missed tip)")
ACN = [(-2.98974021026, 2.415821643), (0.944876236554, 0.8908654569), (-0.643416747409, -1.07019451211)]
for i, av in enumerate(ACN):
    ax.plot(*av, "D", ms=8, mec="k", mfc="indigo", zorder=5,
            label="acnodes (3)" if i == 0 else None)
ax.plot(0.233073, -0.001790, "o", ms=9, mec="k", mfc="lime", zorder=5, label="fold w=1/2")
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC^2")
ax.set_title("(a) REAL atlas of the nonic wall: 1-root sea (pale), 3-root wedge (green),\nno missed cone (odd chamber) - only the whisker's TIP is missed over R", fontsize=10.5)
ax.legend(loc="lower right", fontsize=9)

# (b) census across chambers
ax = axs[0][1]
cens = {}
for d in (5, 6, 7, 8):
    cens[d] = json.load(open(f"atlas{d}_realcensus.json"))
print("censuses:", cens)
labels, vals, cols = [], [], []
for d in (5, 6, 7, 8):
    tot = sum(cens[d].values())
    for k in sorted(cens[d], key=int):
        labels.append(f"d={d}\nfiber {d+1}\n{k} preim")
        vals.append(100*cens[d][k]/tot)
        cols.append("#c23b22" if int(k) == 0 else ("#1f77b4" if int(k) % 2 else "#2ca02c"))
xpos = np.arange(len(labels))
ax.bar(xpos, vals, color=cols)
ax.set_xticks(xpos); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("% of 200k real targets")
ax.set_title("(b) real-census rhyme: even fibers carry the missed CONE (red, ~8.7%),\nodd fibers are R-surjective a.e. ({1,3} split ~ 82.6/17.4 both nights)", fontsize=10.5)

# (c) monodromy
ax = axs[1][0]; ax.axis("off")
mono = open("atlas8_monodromy.txt").read().strip().split("\n")
head, loops = mono[0], mono[1:]
txt = "MONODROMY CERTIFICATE - fiber cover has 9 sheets\n" + "="*52 + "\n"
for ln in loops:
    name, perm = ln.split(": ", 1)
    txt += f"  {name:<34s} {perm}\n"
txt += "\n  generated group:  |G| = 362880  = 9!\n"
txt += "  transpositions in G: 36/36   ->  G = S9  (5th chamber)\n"
txt += "  loop min|D8| kept > 1e15: no wall crossings;\n  4000/8000-step refinements agree on every loop."
ax.text(0.02, 0.96, txt, transform=ax.transAxes, fontsize=10.5, va="top", family="monospace")
ax.set_title("(c) braid of the 9 sheets", fontsize=11)

# (d) term law + fingerprints
ax = axs[1][1]
ns = np.arange(3, 11)
cone = ns*(ns+1)/2 + 1
terms_obs = [5, 9, 14, 20, 26, 34, 43]
holes = np.array([2, 2, 1, 2, 3, 3, 3])
ax.plot(ns, cone, "k--", lw=1.2, label="support cone C(n) = n(n+1)/2 + 1  (THM)")
ax.plot(np.arange(3, 10), terms_obs, "o-", color="crimson", ms=8, label="observed terms (d=2..8)")
ax.plot(np.arange(3, 10), cone[:7] - holes, "x", color="gray", ms=7, mew=2,
        label="cone - holes: constants & s^1 always, r^{n-2} for n>=7")
ax.plot([10], [53], "*", ms=20, mec="k", mfc="gold", zorder=5, label="n=10 prediction: 53")
for n_, t_ in zip(np.arange(3, 10), terms_obs):
    ax.annotate(str(t_), (n_, t_), textcoords="offset points", xytext=(4, -12), fontsize=9)
ax.set_xticks(ns); ax.set_xlabel("fiber degree n")
ax.set_title("(d) wall term law: terms = C(n) - 2 - [n>=7];\ns^9 = -2^58*3^5 fingerprint; magnitude law EXACT in all 7 chambers", fontsize=10.5)
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig("atlas8_atlas.png", dpi=170, bbox_inches="tight")
print("saved atlas8_atlas.png")
