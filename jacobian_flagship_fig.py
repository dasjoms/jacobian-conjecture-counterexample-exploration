"""
NOTE 22 figure: THE THRESHOLD DECK, four panels.
(a) d=2 wall (sextic D2 in (s,r)) with its cusp (1/3, 1/27): the missed cube
(b) d=3 wall D3 with node (-1,-1) and 2 real cusps: the missed square
(c) partition ledger n=3..8: realized / killed-by-X / rescued
(d) THE DECK: chamber cards n=3..13 with verdicts
"""
import sympy as sp, numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

w, s, r = sp.symbols("w s r")
R = sp.Rational

def tower_p(d):
    m = d * (d + 1)
    return sp.expand(2*w - 3*w**2 + w*(1 - w)*(w**(d - 2) - R(6, m)))

def Phi(d):
    return sp.expand(sp.integrate(tower_p(d), w))

# ---- exact walls ----
D2 = sp.expand(sp.discriminant(Phi(2) - s*w + r, w))
h3 = Phi(3) - s*w + r
D3 = sp.factor(sp.resultant(h3, sp.diff(h3, w), w))
print("D2 =", D2)
print("D3 =", D3)

fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 2, left=0.05, right=0.98, top=0.90, bottom=0.06, wspace=0.20, hspace=0.34)

# ================================================================ (a)
ax = fig.add_subplot(gs[0, 0])
xx = np.linspace(-0.7, 1.2, 900); yy = np.linspace(-0.45, 0.6, 900)
S, Rg = np.meshgrid(xx, yy)
D2n = sp.lambdify((s, r), D2, "numpy")
Z = D2n(S, Rg)
ax.contour(S, Rg, Z, levels=[0], colors=["#1a56c4"], linewidths=[1.6])
ax.plot(R(1,3), R(1,27), "*", ms=20, mfc="#ffd21a", mec="k", mew=1.2, zorder=6)
ax.annotate("cusp $(s^*,r^*)=(1/3,\\,1/27)$\nfiber $h=-(w-\\frac{1}{3})^3$ perfect cube\n$\\Rightarrow$ targets $\{(1/(27C^2),\\,1/(3C),\\,C)\\}$ MISSED",
            xy=(1/3, 1/27), xytext=(0.52, 0.34), fontsize=9.5,
            arrowprops=dict(arrowstyle="->", color="k"), ha="left",
            bbox=dict(fc="#fff7d6", ec="gray"))
ax.set_xlabel("$s = BC$", fontsize=12); ax.set_ylabel("$r = AC^2$", fontsize=12)
ax.set_title("(a) chamber $n=3$: wall $D_2 = 4s^3-s^2-18sr+27r^2+4r=0$ (note 1's cubic-dual)\nunique singularity = the cusp; everything else attained", fontsize=11)
ax.grid(alpha=0.25)

# ================================================================ (b)
axb = fig.add_subplot(gs[0, 1])
xx = np.linspace(-3.6, 2.4, 1000); yy = np.linspace(-3.2, 2.6, 1000)
S2, R2 = np.meshgrid(xx, yy)
D3n = sp.lambdify((s, r), D3, "numpy")
m3 = 12
p3 = tower_p(3)
Z2 = D3n(S2, R2)
Z2c = np.sign(Z2) * np.log1p(np.abs(Z2))
axb.contour(S2, R2, Z2c, levels=[0], colors=["#1a56c4"], linewidths=[1.6])
axb.plot(-1, -1, "*", ms=22, mfc="#ffd21a", mec="k", mew=1.2, zorder=6)
axb.annotate("node $(s^*,r^*)=(-1,-1)$, contacts $w\\in\\{1,-2\\}$\nfiber $h=-\\frac{1}{4}(w^2+w-2)^2$ perfect square\n$\\Rightarrow$ ALPÖGE'S MAP misses $\{(-1/C^2,\\,-1/C,\\,C)\\}$",
             xy=(-1, -1), xytext=(-3.4, 1.2), fontsize=9.5,
             arrowprops=dict(arrowstyle="->", color="k"), ha="left",
             bbox=dict(fc="#fff7d6", ec="gray"))
taus = [(-1 - 3**0.5)/2, (-1 + 3**0.5)/2]
for tau in taus:
    sc = float(p3.subs(w, tau)); rc = float((w*p3 - Phi(3)).subs(w, tau))
    axb.plot(sc, rc, "s", ms=9, mfc="none", mec="#c41a3a", mew=2, zorder=5)
axb.annotate("2 real cusps (rescued:\nfibers have simple roots)", xy=(-2.3, 0.85), xytext=(0.3, -1.62),
             fontsize=9, color="#c41a3a", arrowprops=dict(arrowstyle="->", color="#c41a3a"))
axb.set_xlabel("$s = BC$", fontsize=12); axb.set_ylabel("$r = AC^2$", fontsize=12)
axb.set_title("(b) chamber $n=4$: wall $D_3$ (from disc of $\\Phi_3-sw+r$)\nunique node = the ONLY missed $(s,r)$ of the tower", fontsize=11)
axb.grid(alpha=0.25)

# ================================================================ (c)
axc = fig.add_subplot(gs[1, 0])
part_table = [
    (3, [("(3)", "realized", "#e8452c", "fiber = cube: MISSED")]),
    (4, [("(4)", "sqfree", "#7a1ac4", "killed: simple flexes"), ("(2,2)", "realized", "#e8452c", "fiber = square: MISSED")]),
    (5, [("(5)", "sqfree", "#7a1ac4", "killed: simple flexes"), ("(3,2)", "T2", "#e07b1a", "killed: no node-cusp overlap"), ("(2,2,1)", "rescue", "#0a6e2e", "simple root rescues")]),
    (6, [("(6)", "sqfree", "#7a1ac4", ""), ("(4,2)", "sqfree", "#7a1ac4", ""), ("(3,3)", "T2", "#e07b1a", "2 flexes, 1 line"), ("(2,2,2)", "T1", "#1a56c4", "tritangent"), ("(2,2,1,1)", "rescue", "#0a6e2e", "")]),
    (7, [("(7),(5,2),(4,3)", "sqfree", "#7a1ac4", ""), ("(3,2,2),(3,3,1)", "T2", "#e07b1a", ""), ("(2,2,2,1)", "rescue", "#0a6e2e", "any marker 1 rescues")]),
    (8, [("all with part $\\geq 4$ or two 3s", "sqfree/T2", "#7a1ac4", ""), ("(2,2,2,2)", "T1", "#1a56c4", "quadritangent"), ("rest", "rescue", "#0a6e2e", "")]),
]
axc.set_xlim(0, 10.8); axc.set_ylim(-0.8, len(part_table))
axc.axis("off")
axc.set_title("(c) the partition ledger: who kills the powerful fibers?\n"
              "(fiber all-multiple $\\Leftrightarrow$ some target lost)", fontsize=11, loc="left")
for i, (n, parts) in enumerate(part_table):
    yv = len(part_table) - 1 - i
    axc.text(0.1, yv + 0.25, f"$n={n}$", fontsize=11, fontweight="bold")
    xs = 1.3
    for (pname, verdict, col, note_) in parts:
        wbox = 0.75 + 0.16 * len(pname)
        fb = FancyBboxPatch((xs, yv), wbox, 0.52, boxstyle="round,pad=0.07",
                            fc=col, ec="k", alpha=0.88)
        axc.add_patch(fb)
        axc.text(xs + wbox/2, yv + 0.26, pname, fontsize=9, color="white",
                 ha="center", va="center", fontweight="bold")
        if note_:
            axc.text(xs + wbox/2, yv - 0.17, note_, fontsize=7.2, color="#333",
                     ha="center", va="center")
        xs += wbox + 0.25
lab = [("#e8452c", "REALIZED -> target missed"), ("#7a1ac4", "killed by simple-flexes (gcd(p',p'')=1)"),
       ("#e07b1a", "killed by T2 (gcd(Cof, p')=1)"), ("#1a56c4", "killed by T1 (Cof squarefree)"),
       ("#0a6e2e", "rescued by a simple root")]
for j, (col, txt) in enumerate(lab):
    axc.add_patch(FancyBboxPatch((0.15 + j*1.95, -0.62), 0.24, 0.24, boxstyle="round,pad=0.03", fc=col, ec="k"))
    axc.text(0.48 + j*1.95, -0.70, txt, fontsize=6.8, va="center")

# ================================================================ (d)
axd = fig.add_subplot(gs[1, 1])
axd.set_xlim(0, 12.4); axd.set_ylim(-1.1, 6.2)
axd.axis("off")
axd.set_title("(d) THE DECK: chambers n = 3..13(+) drawn as cards", fontsize=11, loc="left")
cards = [
    (3, "#e8452c", ["C-surjective: NO", "misses cube-curve", "wall: 2 cusps", "0 nodes", "R: NO"]),
    (4, "#e8452c", ["C-surjective: NO", "misses sq-curve", "ALPOGE's chamber", "1 node (that one)", "R: NO"]),
    (5, "#0a6e2e", ["C-surjective: YES", "gcd certs T1^T2", "+ simple flexes", "3 cusps, 3 nodes", "R: NO (un-rescue)"]),
    (6, "#0a6e2e", ["C-surjective: YES", "4 cusps, 6 nodes", "R: NO (ghost)"]),
    (7, "#0a6e2e", ["C-surjective: YES", "5 cusps, 10 nodes", "monodromy S_7", "R: NO"]),
    (8, "#0a6e2e", ["C-surjective: YES", "6 cusps, 15 nodes", "15=1+2+12 ok", "R: NO"]),
    (9, "#0a6e2e", ["C-surjective: YES", "7 cusps, 21 nodes", "acnode stair ok", "R: NO"]),
    (10, "#0a6e2e", ["C-surjective: YES", "8 cusps, 28 nodes", "Gal S_10 exact", "R: NO"]),
    (11, "#0a6e2e", ["C-surjective: YES", "9 cusps, 36 nodes", "corner locked", "R: NO"]),
    (12, "#0a6e2e", ["C-surjective: YES", "10 cusps, 45 nodes", "LAST CHAMBER", "R: NO"]),
    (13, "#1a56c4", ["C-surjective: YES*", "11 cusps, 55 nodes", "*corridor only", "(no descent)"]),
    (14, "#999999", ["CONJECTURE: YES", "T1^T2 open here", "the corridor", "remains"]),
]
for i, (n, col, lines) in enumerate(cards):
    x0 = 0.12 + (i % 6) * 2.05
    y0 = 3.0 - (i // 6) * 3.25
    fb = FancyBboxPatch((x0, y0), 1.92, 2.75, boxstyle="round,pad=0.09", fc="white", ec=col, lw=2.2)
    axd.add_patch(fb)
    axd.text(x0 + 0.96, y0 + 2.45, f"n = {n}", fontsize=11, fontweight="bold", ha="center", color=col)
    for k, ln in enumerate(lines):
        axd.text(x0 + 0.13, y0 + 1.95 - 0.5 * k, ln, fontsize=6.6, va="center",
                 fontweight="bold" if ("YES" in ln or "NO" in ln) else "normal",
                 color=col if ("YES" in ln or "NO" in ln) else "#222222")
fig.suptitle("NOTE 22 — THE THRESHOLD DECK: powerful fibers kill surjectivity exactly in chambers $n=3,4$;\n"
             "from $n=5$ the budget law's gcd certificates rescue every fiber (d = 4..12 certified, all d conjectured)",
             fontsize=12.5, y=0.985)
fig.savefig("flagship_figure.png", dpi=150)
print("figure saved")
