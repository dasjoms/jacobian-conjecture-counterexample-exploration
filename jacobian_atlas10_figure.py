"""NOTE 14 figure: the whisker chamber, 4 panels, 2100x1400."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(21.0, 14.0))
plt.subplots_adjust(hspace=0.32, wspace=0.24, top=0.93, bottom=0.08, left=0.06, right=0.985)
fig.suptitle("NOTE 14 — THE WHISKER CHAMBER: fiber 11, wall degree 11, 64 terms, 9 cusps (1 real), "
             "36 nodes (0 crunodes, 4 acnodes), monodromy $S_{11}$", fontsize=15.5, fontweight="bold")

def p10(w): return -w**10 + w**9 - 162*w**2/55 + 107*w/55
def Phi11(w): return -w**11/11 + w**10/10 - 54*w**3/55 + 107*w**2/110
def tau(w): return w*p10(w) - Phi11(w)

# ---------------- (a) whisker atlas ----------------
ax = axes[0][0]
d = np.load("atlas10_grid.npz")
S, R, grid = d["S"], d["R"], d["grid"]
ax.pcolormesh(S, R, grid, cmap="plasma", shading="auto", alpha=0.85)
T = np.linspace(-1.6, 2.8, 3000)
ax.plot([p10(t) for t in T], [tau(t) for t in T], "-", color="#ff1744", lw=2.4, alpha=0.95, label="wall (param)", zorder=6)
rc = (0.321272, 0.035372)
ax.plot(*rc, "*", color="#00ffcc", ms=18, mec="k", mew=0.8, label="1 real cusp (whisker tip)")
acs = [(-3.5335319, 2.4380180), (-0.2244506, 1.7112555), (1.0264806, 0.1411494), (-0.84220146, -1.1018862)]
for i, a in enumerate(acs):
    ax.plot(*a, "D", color="#7bffd4", ms=9, mec="k", mew=0.8,
            label="4 acnodes" if i == 0 else None)
ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) the whisker chamber: real pre-image counts {1,3}; missed set = wall itself")
ax.legend(fontsize=10, loc="lower right")

# ---------------- (b) term law + support map ----------------
ax = axes[0][1]
ns = np.arange(3, 12); terms = [5,9,14,20,26,34,43,53,64]
lawl = ns*(ns+1)/2 - 2
ax.plot(ns, terms, "o", color="#7b3294", ms=10, mec="k", label="measured terms (9 chambers)")
ax.plot(ns, lawl, "-", color="#b8860b", lw=2, label="$n(n+1)/2 - 2$  (three holes, $n \\geq 7$)")
ax.plot([11], [64], "*", color="#d01c1c", ms=17, mec="k", label="64 — locked in note 12/13, tonight GREEN")
ax.set_xlabel("fiber degree n"); ax.set_ylabel("wall monomial count")
ax.legend(fontsize=9.5, loc="upper left"); ax.grid(alpha=0.3)
ax.set_title("(b) TERM LAW: support = cone − {0, $s$, $r^{n-2}$, $r^n$(fict.)} + $s^n$ — EXACT on stored walls")
# inset support map for n=11
axin = ax.inset_axes([0.50, 0.40, 0.46, 0.56])
pts = [(i,j) for i in range(12) for j in range(12) if 10*j+11*i <= 110]
holes = [(0,0), (1,0), (0,9), (0,11)]
keep = [p for p in pts if p not in holes]
axin.scatter([p[1] for p in keep], [p[0] for p in keep], s=26, color="#7b3294", label="present")
axin.scatter([p[1] for p in holes], [p[0] for p in holes], s=70, facecolor="none", edgecolor="#d01c1c", lw=1.6, label="hole")
axin.scatter([0], [11], s=90, marker="s", color="#b8860b", zorder=5, label="$s^{11}$ (out of cone)")
lbl = {(0,0): ("const", (7,-2)), (1,0): ("$s$", (7,-2)), (0,9): ("$r^9$", (-26,2)), (0,11): ("$r^{11}$ fict.", (-50,-10))}
for hh, (nm, off) in lbl.items():
    axin.annotate(nm, (hh[1], hh[0]), textcoords="offset points", xytext=off, fontsize=9, color="#d01c1c", fontweight="bold")
axin.annotate("$s^{11}$ (out of cone)", (0, 11), textcoords="offset points", xytext=(4,-3), fontsize=9, color="#b8860b", fontweight="bold")
axin.set_xlim(-0.7, 12); axin.set_ylim(-0.9, 12.2)
axin.set_xlabel("j (in $r^j$)", fontsize=8); axin.set_ylabel("i (in $s^i$)", fontsize=8)
axin.tick_params(labelsize=7); axin.set_title("support map, $n=11$", fontsize=9)

# ---------------- (c) census + certificate table ----------------
ax = axes[1][0]
cen = json.load(open("atlas10_realcensus.json"))
b1 = [cen.get("1",0)/2000, cen.get("3",0)/2000, cen.get("5",0)/2000]
x0 = np.arange(3)
ax.bar(x0-0.18, b1, width=0.36, color="#7b3294", label="d=10 census (200k): {1:82.65, 3:17.35, 5:0}")
c9 = json.load(open("atlas9_realcensus.json")); tot9 = sum(c9.values())
b9 = [c9.get("0",0)/tot9*100, c9.get("2",0)/tot9*100, c9.get("4",0)/tot9*100]
ax.bar(x0+0.18, b9, width=0.36, color="#b8860b", label="d=9 (cone): {0:8.73, 2:84.03, 4:7.25}")
ax.set_xticks(x0); ax.set_xticklabels(["missed-0/1\n(whisker vs cone)", "missed-2/3", "missed-4/5"])
ax.set_ylabel("%"); ax.legend(fontsize=9)
ax.set_title("(c) parity census: odd fibers never cone — the missed set is measure-zero (whisker)")
# table of the eliminant arithmetic
txt = ("eliminant deg 90 = (55·p₁₀′)² · cofactor[72]   K = 55² = 3025 EXACT\n"
       "cofactor squarefree + coprime to p₁₀′ (exact gcds over ℚ)\n"
       "budget (n−1)(n−2) = 90 = 2·9 + 2·36 ✓   (9 cusps, 36 nodes)\n"
       "nodes: 0 crunodes (dance) + 4 acnodes (staircase) + 32 complex\n"
       "acnode residuals: 1/7 real roots each (odd septic)\n\n"
       "fingerprint [s¹¹]:  −2²⁰·5²⁸·11⁹   sign = (−1)ⁿ  (9/9)\n"
       "content(11) = 2¹¹·5³·11²   support ⊆ primes(n(n−1)) = {2,5,11}\n"
       "magnitude law ratio = 1   (L = 110)\n"
       "gcd(p₁₀′, p₁₀″) = 1   → 9 simple cusps;  Sturm: 1 real (dance 1,2,…,1)")
ax.text(0.02, -0.62, txt, transform=ax.transAxes, fontsize=11, family="monospace",
        bbox=dict(boxstyle="round", fc="#f5f5ff", ec="#7b3294"))

# ---------------- (d) monodromy + scoreboard ----------------
ax = axes[1][1]; ax.axis("off")
mono = open("atlas10_monodromy.txt").read().splitlines()
lines = ["MONODROMY — Jordan certificate for $S_{11}$:"] + [f"   {L}" for L in mono[1:]] + [
    "", "   transitive: True; 2-homogeneous: 55/55 pairs; transposition (9,10) present",
    "   => G = S11 (Jordan).   bonus |G| ≥ 400006    ── scoreboard: 18/18 GREEN ──"]
sb = [
    ("terms(11) = 64; holes {const, s, r⁹}", "GREEN"),
    ("D10 irreducible, deg 11, param identity", "GREEN"),
    ("s² | D10(s,0); D10(0,0) = 0", "GREEN"),
    ("magnitude law L=110, ratio 1", "GREEN"),
    ("content = 2¹¹·5³·11²; support ⊆ {2,5,11}", "GREEN"),
    ("fingerprint sign (−1)¹¹, 9/9 chamber run", "GREEN"),
    ("9 cusps; Sturm 1 real; NOT HIT (whisker)", "GREEN"),
    ("eliminant = (55p′)²·cofactor[72], K=3025", "GREEN"),
    ("cofactor squarefree + coprime", "GREEN"),
    ("36 nodes = 0 crunodes + 4 acnodes + 32 cplx", "GREEN"),
    ("budget 90 = 18 + 72", "GREEN"),
    ("fibers generic/fold/cusp/acnode 11/9/8/7", "GREEN"),
    ("slopes fold 1/2, cusp 2/3", "GREEN"),
    ("census {1:[80.5,84.5],3:[15.5,19.5],5:[0,0.5]}", "GREEN 82.65/17.35/0"),
    ("hinge antipodal (p₁=2q₂ identity, K=d)", "GREEN ±0.90603"),
    ("det JF = 1 at 5/5 rational points", "GREEN"),
    ("monodromy S₁₁ (Jordan)", "GREEN"),
    ("refined term law exact on stored walls n=10,11", "GREEN"),
]
tbl = ax.table(cellText=[["#", "prediction", "verdict"]] + [[str(i+1), a, b] for i, (a, b) in enumerate(sb)],
               colWidths=[0.06, 0.64, 0.30], loc="lower center", cellLoc="left")
tbl.auto_set_font_size(False); tbl.set_fontsize(7.8); tbl.scale(1, 1.05)
for j in range(3): tbl[0, j].set_facecolor("#333333"); tbl[0, j].set_text_props(color="w", weight="bold")
for i in range(1, len(sb)+1): tbl[i, 2].set_facecolor("#e6ffe6")
ax.text(0.0, 0.995, "\n".join(lines), transform=ax.transAxes, fontsize=10.5, va="top", family="monospace")
ax.set_title("(d) 18/18 green — the eleventh chamber obeys every law at once")
fig.savefig("atlas10_figure.png", dpi=100)
print("saved atlas10_figure.png")
