"""Note 20 figure: the last chamber atlas (2x2 gridspec, manual spacing)."""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

d = 11; m = d*(d+1)
def p11v(w): return -w**d + w**(d-1) - (3-6/m)*w**2 + (2-6/m)*w
def Phi11v(w): return -w**(d+1)/(d+1) + w**d/d - (3-6/m)*w**3/3 + (2-6/m)*w**2/2

fig = plt.figure(figsize=(11.5, 9.2), dpi=150)
gs = fig.add_gridspec(2, 2, left=0.07, right=0.975, top=0.93, bottom=0.075, wspace=0.24, hspace=0.34)

# ---------- (a) the d=11 real atlas ----------
axa = fig.add_subplot(gs[0,0])
smn, smx, rmn, rmx = -3.7, 1.6, -1.35, 2.7
G = 70
S0 = np.linspace(smn, smx, G); R0 = np.linspace(rmn, rmx, G)
co = np.zeros(d+2); co[0] = -1/(d+1); co[1] = 1/d; co[d-2] = -(1-2/m); co[d-1] = (1-3/m)
pc = np.array([-1.0, 1.0, 0,0,0,0,0,0,0, -(3-6/m), (2-6/m)])
cnt = np.zeros((G, G), int)
for i, sv in enumerate(S0):
    for j, rv in enumerate(R0):
        cc = co.copy(); cc[-2] = -sv; cc[-1] = rv
        rt = np.roots(cc)
        cnt[j, i] = sum(1 for z in rt if abs(z.imag) < 1e-8 and abs(sv - np.polyval(pc, z.real)) > 1e-9)
cmap = ListedColormap(["#f6d9d9", "#fdf6e3", "#dcebf7", "#cfe3d3", "#cfe3d3", "#cfe3d3", "#cfe3d3"])
axa.imshow(np.clip(cnt, 0, 6), origin="lower", extent=[smn, smx, rmn, rmx], cmap=cmap,
           aspect="auto", alpha=0.75, zorder=0)
t = np.linspace(-3.35, 3.35, 6000)
axa.plot(p11v(t), t*p11v(t) - Phi11v(t), color="#8b0000", lw=1.6, zorder=3, label="wall $t\\mapsto(s,r)$")
marks = [
    (( 0.3232622472,  0.03564372163), "*", "#d4a500", 130, "cusps", 1),
    ((-3.491074673,   1.684483226),   "*", "#d4a500", 130, None, 1),
    ((-0.9416193117, -0.9419704421),  "s", "#c2185b", 75,  "corner (crunode)", 1),
    ((-2.468169393,   2.386341511),   "D", "#0288d1", 46,  "acnodes", 1),
    (( 0.5466482194,  1.2868017516),  "D", "#0288d1", 46,  None, 1),
    (( 0.837322656,  -0.0881120911),  "D", "#0288d1", 46,  None, 1),
    ((-0.8937318943, -1.1033381706),  "D", "#0288d1", 46,  None, 1),
    ((-1.0, -1.0), "P", "#1b5e20", 85, "pin image $(-1,-1)$", 1),
]
import mpmath as mp
mp.mp.dps = 50
w_ = mp.mpc(1)
tstar = mp.mpf("-1.0834586214993059682")
sstar = -tstar**11 + tstar**10 - mp.mpf(65)/22*tstar**2 + mp.mpf(43)/22*tstar
phistar = -tstar**12/12 + tstar**11/11 - mp.mpf(65)/66*tstar**3 + mp.mpf(43)/44*tstar**2
rstar = tstar*sstar - phistar
marks.append(((float(sstar), float(rstar)), "h", "#00695c", 64, "shadow image", 1))
seen_lab = set()
for (xx, yy), mk, cl, sz, lab, z in marks:
    axa.scatter([xx], [yy], marker=mk, s=sz, c=cl, zorder=5 if mk != "P" else 6,
                edgecolors="k", linewidths=0.5, label=lab if lab and lab not in seen_lab else None)
    if lab: seen_lab.add(lab)
axa.set_xlim(smn, smx); axa.set_ylim(rmn, rmx)
axa.set_xlabel("$s$"); axa.set_ylabel("$r$")
axa.set_title("(a)  chamber $n=12$ real atlas: 0/2/4-real regions, wall, entourage", fontsize=10)
axa.legend(loc="lower right", fontsize=7.5, framealpha=0.95)

# ---------- (b) support cone ----------
axb = fig.add_subplot(gs[0,1])
NN = 12
holes = {(0,0), (1,0), (0,10)}
inc = []
for i in range(NN+1):
    for j in range(NN+1):
        if (NN-1)*i + NN*j <= NN*(NN-1):
            inc.append((i, j))
for (i, j) in inc:
    c = "#4caf50" if (i, j) not in holes else "#e53935"
    ec = "k" if (i, j) in holes else "none"
    axb.add_patch(plt.Rectangle((i-0.45, j-0.45), 0.9, 0.9, facecolor=c, edgecolor=ec, lw=1.1,
                                alpha=0.55 if (i,j) not in holes else 0.9))
xx = np.linspace(0, 12, 50)
axb.plot(xx, (NN*(NN-1) - (NN-1)*xx)/NN, "k--", lw=1.1)
for (i,j), lab, off in [((0,0),"const",(8,8)), ((1,0),"$s^1$",(8,-14)), ((0,10),"$r^{10}$",(10,-12))]:
    axb.annotate(lab, (i, j), textcoords="offset points", xytext=off, fontsize=9, color="#b71c1c")
axb.annotate("$s^{12}$", (12, 0), textcoords="offset points", xytext=(-32, 10), fontsize=9, color="#1b5e20")
axb.annotate("$r^{11}$", (0, 11), textcoords="offset points", xytext=(10, 2), fontsize=9, color="#1b5e20")
axb.set_xlim(-1, 13); axb.set_ylim(-1, 12.4)
axb.set_xlabel("$i$ (power of $s$)"); axb.set_ylabel("$j$ (power of $r$)")
axb.set_title("(b)  support of $D_{11}$: 76 green seats, 3 red holes, cone $11i+12j\\leq132$", fontsize=10)

# ---------- (c) census migration ----------
axc = fig.add_subplot(gs[1,0])
fresh = json.load(open("atlas11_stageD.json"))["fresh_census_pct"]
dd_even = [4,6,8,10]; row3 = [fresh[str(k)].get("3",0.0) if "3" in fresh[str(k)] else fresh[str(k)].get(3,0.0) for k in dd_even]
dd_odd  = [5,7,9,11]; row0 = [fresh[str(k)]["0"] if "0" in fresh[str(k)] else fresh[str(k)].get(0,0.0) for k in dd_odd]
row4    = [fresh[str(k)]["4"] if "4" in fresh[str(k)] else fresh[str(k)].get(4,0.0) for k in dd_odd]
row3 = [fresh[str(k)].get("3", fresh[str(k)].get(3,0.0)) for k in dd_even]
row0 = [fresh[str(k)].get("0", fresh[str(k)].get(0,0.0)) for k in dd_odd]
row4 = [fresh[str(k)].get("4", fresh[str(k)].get(4,0.0)) for k in dd_odd]
axc.plot(dd_even, row3, "o-", color="#6a1b9a", ms=6, lw=1.5, label="3-real row (whisker: $n$ odd)")
axc.plot(dd_odd, row0, "s-", color="#00695c", ms=6, lw=1.5, label="cone 0-real row ($n$ even)")
axc.plot(dd_odd, row4, "^-", color="#ad1457", ms=6, lw=1.5, label="4-real row ($n$ even)")
axc.axhspan(8.67, 8.84, xmin=0.86, xmax=0.985, color="#00695c", alpha=0.18)
axc.plot([11.0], [8.7556], marker="_", ms=22, mew=2.4, color="#004d40")
axc.annotate("published lock [8.67, 8.84];\nlaw center 8.7556", (11, 8.7556), textcoords="offset points",
             xytext=(-170, 26), fontsize=8, color="#004d40")
lab_off = {5: (-12, -16), 7: (-12, -16), 9: (-12, -16), 11: (-52, -14)}
for k, v in zip(dd_even, row3): axc.annotate(f"{v:.3f}", (k, v), textcoords="offset points", xytext=(-6, 9), fontsize=7.5, color="#6a1b9a")
for k, v in zip(dd_odd, row0): axc.annotate(f"{v:.3f}", (k, v), textcoords="offset points", xytext=lab_off[k], fontsize=7.5, color="#00695c")
for k, v in zip(dd_odd, row4): axc.annotate(f"{v:.3f}", (k, v), textcoords="offset points", xytext=(9, 2), fontsize=7.5, color="#ad1457")
axc.set_xlim(3.4, 12.1); axc.set_ylim(5.6, 18.2)
axc.set_xlabel("$d = n-1$"); axc.set_ylabel("census %  (240k, $\\sigma=1.5$)")
axc.set_title("(c)  parity census migration across the tour", fontsize=10)
axc.legend(loc="center left", fontsize=8, bbox_to_anchor=(0.02, 0.62))

# ---------- (d) corner chases shadow ----------
axd = fig.add_subplot(gs[1,1])
sh = [(3,-1.0), (5,-0.8825387208), (7,-0.9098546437), (9,-0.9289658123), (11,-0.9417878440), (13,-0.9508166086)]
cr = [(5,-0.88188341), (7,-0.9094571137), (9,-0.9287162891), (11,-0.9416193117)]
crr = [(5,-0.88337174), (7,-0.9103141728), (9,-0.9292425119), (11,-0.9419704421)]
axd.axhline(-1, color="#888", ls=":", lw=1.2)
axd.annotate("the pin: $-1$", (6.3, -1.0007), fontsize=8.5, color="#555")
axd.plot([q[0] for q in sh], [q[1] for q in sh], "o-", color="#00695c", ms=6, lw=1.6, label="shadow $s^*(d)$")
axd.plot([q[0] for q in cr], [q[1] for q in cr], "s-", color="#c2185b", ms=7, lw=1.6, label="crunode $s_c(d)$")
axd.plot([q[0] for q in crr], [q[1] for q in crr], "^--", color="#ef6c00", ms=6, lw=1.3, label="crunode $r_c(d)$")
axd.scatter([3], [-1.0], marker="P", s=110, c="#1b5e20", zorder=6, edgecolors="k", linewidths=0.5)
axd.annotate("$d=3$ coalescence:\ncorner $=$ shadow $=$ pin", (3, -1.0), textcoords="offset points", xytext=(14, 6), fontsize=8, color="#1b5e20")
gaps = [(5, 6.5531e-4), (7, 3.9753e-4), (9, 2.4952e-4), (11, 1.6853e-4)]
axi = axd.inset_axes([0.52, 0.575, 0.45, 0.38])
axi.loglog([q[0] for q in gaps], [q[1] for q in gaps], "o", color="#4527a0", ms=5)
xx = np.linspace(4.8, 11.5, 60)
axi.loglog(xx, 0.0204/xx**2, "-", color="#4527a0", lw=1.2, alpha=0.7, label="$0.0204/d^2$")
axi.set_title("gap $s_c-s^*$", fontsize=8)
axi.tick_params(labelsize=7); axi.legend(fontsize=7, loc="lower left", framealpha=0.9)
axi.grid(True, which="both", lw=0.3, alpha=0.4)
for k, v, dx, dy in [(5, -0.88188341, 8, 11), (7, -0.9094571137, -58, 4), (9, -0.9287162891, 8, -16), (11, -0.9416193117, 8, -16)]:
    axd.annotate(f"{v:.5f}", (k, v), textcoords="offset points", xytext=(dx, dy), fontsize=7.3, color="#880e4f")
axd.set_xlim(2.6, 13.4); axd.set_ylim(-1.006, -0.845)
axd.set_xlabel("$d$"); axd.set_ylabel("$s$")
axd.set_title("(d)  the corner chases the shadow (both $\\to$ the pin)", fontsize=10)
axd.legend(loc="lower right", fontsize=7.8)

fig.suptitle("chamber $n=12$ ($d=11$): 76 terms, $K=484$, 10 cusps (2 real, both hit), 45 nodes $=1\\oplus 4\\oplus 40$," +
             " corner $s_c=-0.94161931$ — the last climb", fontsize=11.5)
fig.savefig("atlas11_figure.png")
print("saved atlas11_figure.png")
