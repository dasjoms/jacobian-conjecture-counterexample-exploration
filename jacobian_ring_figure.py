"""NOTE 13 figure: the ringing constant, 4 panels on one 2100x1400 canvas."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import sympy as sp

mp.mp.dps = 50
sig = 1.5
mstar = 0.08644464653779008
A0, B0 = 0.0311269250795446, 0.0185094013627119

a2 = json.load(open("ring_stageA2.json")); sc = json.load(open("ring_stageC.json"))
sweep = {int(k): v*100 for k, v in a2["sweep"].items()}
exact_full = {int(k): v*100 for k, v in sc["exact_full"].items()}
fit = sc["fit"]; Aeff, Bfit = fit["Aeff"], fit["Bfit"]

fig, axes = plt.subplots(2, 2, figsize=(21.0, 14.0))
plt.subplots_adjust(hspace=0.30, wspace=0.22, top=0.93, bottom=0.07, left=0.06, right=0.985)
fig.suptitle("NOTE 13 — THE RINGING CONSTANT:  $m(d) = m^* + [A - B(\\ln d - \\ln\\ln d)]/d$,  "
             "with $m^* = 8.64446465377901\\%$", fontsize=17, fontweight="bold")

T = [int(v) for v in sorted(exact_full)]; TM = [exact_full[d] for d in T]
S_ = sorted(sweep); SM = [sweep[d] for d in S_]
lawd = np.linspace(40, 220, 800)
Lw_ = np.log(lawd) - np.log(np.log(lawd))
lawM = 100*(mstar + (Aeff - Bfit*Lw_)/lawd)
lawM0 = 100*(mstar + (A0 - B0*Lw_)/lawd)

# ---------------- (a) the ring ----------------
ax = axes[0][0]
ax.axhline(100*mstar, color="k", lw=1.3, ls="-", label="$m^*$ = 8.644465% (exact integral)")
ax.plot(S_, SM, "-", color="#888888", lw=1.4, alpha=0.8, label="GOLD model sweep (dense)")
ax.plot(lawd, lawM, "-", color="#b8860b", lw=2.2, label=f"fitted law: A={Aeff:.4f}, B={Bfit:.4f}")
ax.plot(T, TM, "o", color="#7b3294", ms=8, mec="k", mew=0.7, label="exact envelope masses (20 chambers)")
ax.plot([161], [exact_full[161]], "*", color="#d01c1c", ms=17, mec="k", mew=0.8,
        label=f"d=161: rebound +2.15e-6 vs d=121  (locked in note 12)")
ax.annotate("overshoot\n8.7736", xy=(9, exact_full[9]), xytext=(16, 8.756), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("crossing $d_0 \\approx 46$", xy=(46, 100*mstar), xytext=(30, 8.632), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("minimum  $d^* \\approx 124\\!-\\!128$\n8.63254", xy=(126, 8.63254), xytext=(150, 8.6375), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k"))
ax.annotate("rebound:\n$m(171) > m(161)$, $m(201) > m(161)$", xy=(196, 8.6342), xytext=(168, 8.6300), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k"))
ax.set_xlabel("d  (tower level; fiber degree n = d+1, cone chambers have odd d)")
ax.set_ylabel("missed-cone mass  [%]")
ax.set_title("(a) the tower rings home: overshoot -> undershoot -> rebound -> $m^*$")
ax.legend(fontsize=9.5, loc="upper right"); ax.grid(alpha=0.3)
ax.set_ylim(8.628, 8.78)

# ---------------- (b) law extraction ----------------
ax = axes[0][1]
dE = np.array([d for d in T if d >= 61], float)
Lk = np.log(dE) - np.log(np.log(dE))
Yk = [(exact_full[int(d)]/100 - mstar)*d for d in dE]
ax.plot(Lk, Yk, "o", color="#7b3294", ms=9, mec="k", label="exact: $d\\,(m(d)-m^*)$")
Lg = np.linspace(2.2, 3.9, 50)
ax.plot(Lg, Aeff - Bfit*Lg, "-", color="#b8860b", lw=2.2,
        label=f"fit: {Aeff:.4f} $-$ {Bfit:.4f}$\\,L$   rms $1.2\\times10^{{-5}}$")
ax.plot(Lg, A0 - B0*Lg, "--", color="#008837", lw=2.0,
        label=f"first principles: A0={A0:.4f}, B0={B0:.4f} (two quadratures)")
ax.text(2.32, 0.030, f"$B_{{fit}}/B_0$ = {Bfit/B0:.4f}", fontsize=14,
        bbox=dict(boxstyle="round", fc="#eeffee", ec="#008837"))
ax.set_xlabel("$L(d) = \\ln d - \\ln\\ln d$")
ax.set_ylabel("$d\\,(m(d) - m^*)$")
ax.set_title("(b) the log-correction: predicted slope within 1.4%")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# ---------------- (c) layer battle at d = 201 ----------------
ax = axes[1][0]
# theory curves (light) + measured min-envelope points on a coarse grid
ALf = lambda u: (u+5)*(1 - mp.log((u+5)/2))
BLs = sp.sympify(json.load(open("ring_probe_sym.json"))["BL"])
BLsy = sp.lambdify(sp.symbols("s"), BLs, "mpmath")
import importlib.util
spec = importlib.util.spec_from_file_location("r3", "/home/user/jacobian_ring_3.py")
r3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(r3)
pcn, tcn, coef = r3.poly_data(201)
sg = np.linspace(-3.4, 0.75, 46)
envE = np.array([r3.env_at(float(sv), pcn, tcn, coef) for sv in sg])
envstar = np.minimum(-sg-2, sg)
ax.plot(sg, 201*(envE-envstar), "o", color="#7b3294", ms=5, mec="k", mew=0.4,
        label="measured $d\\cdot(\\mathrm{env}-\\mathrm{env}^*)$, d = 201")
# model pieces
def c_of(u, d):
    lam = float(d*abs(u + 1))          # |s|-1 up the sign: c > 0
    if lam < 1e-6: return lam/5
    a, b = 0.0, max(2.0, math.log(lam) + 2.0)   # bracket: c(4+e^c) monotone
    for _ in range(80):
        m = (a+b)/2
        if m*(4+math.exp(m)) < lam: a = m
        else: b = m
    return (a+b)/2
sR = np.linspace(-3.4, -1.0, 200)
dmR = [201*((c_of(sv, 201)*(sv+1))/201 + (c_of(sv, 201)-1)*math.exp(c_of(sv, 201))/201**2) for sv in sR]
sL = np.linspace(-3.4, 0.75, 400)
dmL = [float(ALf(mp.mpf(float(sv)))) + float(BLsy(mp.mpf(float(sv))))/201 for sv in sL]
ax.plot(sL, dmL, "-", color="#008837", lw=2.0, label="left layer: $A_L(s) + B_L(s)/d$")
ax.plot(sR, dmR, "-", color="#b8860b", lw=2.0, label="right layer: $c(s{+}1) + (c{-}1)e^c/d$")
ax.axvline(2*math.e - 5, color="#d01c1c", lw=1.2, ls=":")
ax.text(-3.3, 1.6, "$s = 2e-5 = 0.4366$: $A_L$ zero\n(the branch later crosses its own line)", fontsize=10, color="#d01c1c")
ax.axhline(0, color="k", lw=0.8, alpha=0.4)
ax.set_xlabel("s"); ax.set_ylabel("$d\\,\\delta\\mathrm{env}$"); ax.set_ylim(-10.5, 2.6)
ax.set_title("(c) the two scales made visible: left $A_L/d \\geq 0$ vs right $-c(|s|{-}1)/d < 0$")
ax.legend(fontsize=9.5, loc="upper left"); ax.grid(alpha=0.3)

# ---------------- (d) scoreboard ----------------
ax = axes[1][1]; ax.axis("off")
rows = [
 ("Prediction (locked before stage-B measurement)", "verdict"),
 ("seven exact masses land in model windows (131..201)", "7/7 GREEN"),
 ("note-12 lock:  m(161) > m(121) = 8.63289%", "+2.150e-6 GREEN"),
 ("slope locks:  m(171) > m(161)-1.5e-6;  m(201) > m(161)-1.5e-6", "+1.353e-6 / +5.781e-6 GREEN"),
 ("layer law $\\tau = -s-2 + A_L/d + B_L/d^2$, symbolically + 100-digit audits", "certified GREEN"),
 ("$A_L(s) = (s+5)(1-\\ln((s+5)/2))$  zero at $s = 2e-5$", "exact GREEN"),
 ("fit slope vs first-principles $B_0$:  within 5%", "1.38% GREEN"),
 ("law min $d^* \\in (119,131)$", "128.3 by fit / 123.6 by sweep GREEN"),
 ("crossing $d_0 \\in (45, 61)$", "45.7 by fit GREEN"),
 ("model tightness d>=45: <= 0.0008pp", "0.00109pp AMBER (= 2nd-order bias)"),
 ("model tightness d in {29,35}: <= 0.004pp", "0.00115pp GREEN"),
 ("2nd-order bias shape $\\propto (\\ln d)^2/d^2$", "flat 0.25-0.28 ratio GREEN"),
]
tbl = ax.table(cellText=rows, colWidths=[0.62, 0.38], loc="center", cellLoc="left")
tbl.auto_set_font_size(False); tbl.set_fontsize(10.5); tbl.scale(1, 1.55)
for j, txt in enumerate(rows[0]):
    tbl[0, j].set_facecolor("#333333"); tbl[0, j].set_text_props(color="w", weight="bold")
for i in range(1, len(rows)):
    c = "#e6ffe6" if "GREEN" in rows[i][1] else "#fff3d0"
    tbl[i, 1].set_facecolor(c)
ax.set_title("(d) scoreboard: 10/11 green, 1 amber — the amber has a theory")
fig.savefig("ring_figure.png", dpi=100)
print("saved ring_figure.png")
