"""
Note 12 figure: THE LIMIT SHADOW.
(a) envelope convergence: env_d(s) for d=9,29,61 vs the limit model min(s, -s-2),
    with the EXACT pin at (-1,-1) holding for every finite d.
(b) missed-cone mass: m(d) (exact envelope masses, no map MC) descending to
    m* = 8.64446%, censuses overlaid, naive shadow 16.29% marked as refuted;
    rebound prediction annotated.
(c) wall-param convergence: (p_d(t), tau_d(t)), t in [-0.98, 1.0], d = 5,7,9 vs
    the cubic param curve; window [-4.5,4.5]^2.
(d) scoreboard + detective-log strip.
"""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp, mpmath as mp
from sympy import symbols, expand, integrate, Rational as R

w = symbols("w"); mp.mp.dps = 60
def seed(d): return expand(2*w - 3*w**2 + w*(1-w)*(w**(sp.Integer(d)-2) - R(6, d*(d+1))))

def poly_data(d):
    p = seed(d); Phi = expand(integrate(p, w))
    pcn = np.array([float(c) for c in sp.Poly(p, w).all_coeffs()])
    tcn = np.array([float(c) for c in sp.Poly(sp.expand(w*p - Phi), w).all_coeffs()])
    coef = [mp.mpf(str(c)) for c in sp.Poly(p, w).all_coeffs()]
    return pcn, tcn, coef

def env_at(sv, pcn, tcn, coef):
    pc = pcn.copy(); pc[-1] -= sv
    coef_s = list(coef); coef_s[-1] = coef_s[-1] - mp.mpf(float(sv))
    f = lambda z: mp.polyval(coef_s, z)
    df = lambda z: mp.polyval([c*(len(coef_s)-1-i0) for i0, c in enumerate(coef_s[:-1])], z)
    taus = []
    for r0v in np.roots(pc):
        if abs(r0v.imag) > 1e-2: continue
        z = mp.mpc(mp.mpf(float(r0v.real)), mp.mpf(float(r0v.imag)))
        for _ in range(25):
            try: dz = f(z)/df(z)
            except ZeroDivisionError: break
            z -= dz
            if abs(dz) < mp.mpf(10)**-45: break
        if abs(mp.im(z)) < mp.mpf(10)**-10 and abs(mp.re(z)) < 8:
            taus.append(float(np.polyval(tcn, float(mp.re(z)))))
    return min(taus) if taus else np.inf

fig = plt.figure(figsize=(15.5, 10.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.34, wspace=0.24)

# ---------- (a) envelope convergence ----------
ax = fig.add_subplot(gs[0, 0])
S = np.arange(-4, 3.01, 0.15)
for d, col, lw in [(9, "#d32f2f", 2.4), (29, "#fb8c00", 1.9), (61, "#1e88e5", 1.5)]:
    D = poly_data(d)
    E = [env_at(float(sv), *D) for sv in S]
    ax.plot(S, E, color=col, lw=lw, label=f"envelope d={d}", zorder=4)
ax.plot(S[S < -1], S[S < -1], "k:", lw=2, label=r"limit: $r=s$  ($s<-1$)")
ax.plot(S[S >= -1], -S[S >= -1] - 2, "k--", lw=2, label=r"limit: $r=-s-2$  ($s>-1$)")
ax.plot([-1], [-1], "*", ms=22, mfc="gold", mec="k", zorder=8,
        label="corner pin $(-1,-1)$: EXACT at every d")
ax.set_xlim(-4, 3); ax.set_ylim(-5.2, 0.6)
ax.set_xlabel("s"); ax.set_ylabel(r"envelope $r_{\rm env}(s)$")
ax.set_title("(a) missed-cone envelope: convergence to min(s, -s-2), pin EXACT", fontsize=11)
ax.legend(fontsize=8.5, loc="lower left")

# ---------- (b) mass trajectory ----------
ax = fig.add_subplot(gs[0, 1])
fin = json.load(open("limits_final.json"))
ms = {int(k): v for k, v in fin["m"].items()}
ds = sorted(ms); ys = [ms[d] for d in ds]
mstar = float(fin["mstar_pct"])
ax.plot(ds, ys, "o-", color="#1b5e20", ms=4.5, lw=1.8, label="exact envelope mass m(d)  [no map-MC]")
for dc, cen, mk in [(5, 8.6900, "s"), (7, 8.7400, "s"), (9, 8.7260, "s")]:
    ax.plot([dc], [cen], mk, ms=11, mfc="none", mec="#7b1fa2", mew=2.2, zorder=7)
ax.plot([], [], "s", ms=9, mfc="none", mec="#7b1fa2", mew=2.2, label="map censuses (200k)")
ax.axhline(16.2888, color="#b71c1c", ls="--", lw=1.6)
ax.text(13, 15.7, "naive cubic shadow 16.29%  -- REFUTED (rescue branches!)", color="#b71c1c", fontsize=9)
ax.axhline(mstar, color="#1565c0", ls="-", lw=2, alpha=0.75)
ax.text(8, mstar+0.12, f"m* = {mstar:.5f}%  exactly (Phi(-1)^2/2 + integral)", color="#1565c0", fontsize=9.5)
ax.annotate("undershoot + rebound\npredicted (locked: m(161) > m(121))",
            xy=(121, ys[-1]), xytext=(46, 8.98),
            arrowprops=dict(arrowstyle="->", color="#455a64"), fontsize=9, color="#455a64")
ax.set_xscale("log"); ax.set_xlim(4, 300)
ax.set_xticks([5, 10, 20, 40, 60, 121, 300]); ax.set_xticklabels(["5", "10", "20", "40", "60", "121", "300"])
ax.set_ylim(8.45, 17.2)
ax.set_ylabel("missed-cone mass (%)"); ax.set_xlabel("chamber d (odd), log scale")
ax.set_title("(b) the constant explained: m(d) -> m* = 8.64446% (+layer corrections)", fontsize=11)
ax.legend(fontsize=8.5, loc="upper right")

# ---------- (c) wall param convergence ----------
ax = fig.add_subplot(gs[1, 0])
TT = np.linspace(-0.98, 1.0, 5001)
for d, col, lw in [(5, "#ef9a9a", 1.2), (7, "#e57373", 1.6), (9, "#d32f2f", 2.2)]:
    p = seed(d); Phi = expand(integrate(p, w))
    pf = sp.lambdify(w, p, "numpy"); phif = sp.lambdify(w, Phi, "numpy")
    ax.plot(pf(TT), TT*pf(TT)-phif(TT), color=col, lw=lw, label=f"chamber d={d}")
pinf = 2*w-3*w**2; Phiinf = expand(integrate(pinf, w))
pf0 = sp.lambdify(w, pinf, "numpy"); phif0 = sp.lambdify(w, Phiinf, "numpy")
ax.plot(pf0(TT), TT*pf0(TT)-phif0(TT), "k--", lw=2.6, label="CUBIC LIMIT (fiber-3 wall's curve)")
ax.plot([-1], [-1], "*", ms=18, mfc="gold", mec="k", zorder=8)
ax.axvspan(-4.5, 4.5, color="#eceff1", alpha=0.35, zorder=0)
ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.5, 4.5)
ax.set_xlabel("s = p(t)"); ax.set_ylabel(r"r = t p(t) - Phi(t)")
ax.set_title("(c) the tower's wall param curves converge to the ORIGINAL cubic", fontsize=11)
ax.legend(fontsize=8.5, loc="lower right")

# ---------- (d) scoreboard strip ----------
ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
rows = [
    ("THE LIMIT SHADOW - detective scoreboard", "k", True),
    ("envelope-mass machinery reproduces censuses (no map MC): devs 0.10/0.03/0.05%", "g", False),
    ("real wall in window  ->  cubic param curve (compacta |t|<=0.95, ratio 1.32)", "g", False),
    ("left layer x: (-5+2e^x, 3-2e^x) = LINE r=-s-2;  right layer env = s", "g", False),
    ("corner pin (-1,-1) EXACT at every d audited (env(-1) = -1.0000)", "g", False),
    ("m* = 0.086444647  exactly; censuses = m* + layer corrections", "g", False),
    ("partA exact: int phi Phi = Phi^2/2 (12-digit check)", "g", False),
    ("algebra: tau_+ - LINE = -2 (t+1)^2 (t-1)  >= 0  on t <= 1", "g", False),
    ("", "k", False),
    ("refuted en route (with fixes):", "k", True),
    ("naive cubic shadow I = 16.2888% != 8.73% - the rescue binds", "r", False),
    ("Newton polish to wrong poly (p vs p-s): invisible at s=0, caught by audit grid", "r", False),
    ("hand Simpsons with dropped rows (twice) - machine cross-check saved us", "r", False),
    ("", "k", False),
    ("locked for next rounds: m(161) > m(121) (rebound); corner(11) s in (-0.950,-0.937)", "#455a64", False),
    ("census(12th) in [8.63%,8.88%]; content-chain law; t2-migration series", "#455a64", False),
]
y = 0.99
for txt, c, bold in rows:
    ax.text(0.02, y, txt, fontsize=10 if bold else 9.2, fontweight="bold" if bold else "normal",
            color={"g": "#1b5e20", "r": "#b71c1c", "k": "black", "#455a64": "#455a64"}[c],
            va="top", transform=ax.transAxes)
    y -= 0.068
ax.set_title("(d) scoreboard", fontsize=11)

fig.suptitle("JACOBIAN LAB - Note 12: THE LIMIT SHADOW - the tower dreams of the cubic, and 8.64% is a theorem",
             fontsize=13.5, fontweight="bold")
fig.savefig("jacobian_limits_figure.png", dpi=150, bbox_inches="tight")
print("saved jacobian_limits_figure.png")
