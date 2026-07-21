"""
Note 21 figure: THE BUDGET LAW in four panels.
  (a) the incidence chamber at d=5: {eq1=0} blue, {te2=0} red, diagonal
      dashed; I = cusps + ordered nodes; counts annotated
  (b) germ zoom at a real flex tau: three branches, slopes -1, -2, +1
  (c) Gamma_5 self-portrait: graph of Phi_5, 2 real flex tangents,
      crunode bitangent (solid), acnode line (dashed, contacts nonreal)
  (d) evidence column d=4..12: stacked budget bars 2(d-1) + (d-1)(d-2)
      = d(d-1), K = den^2 diamonds (log scale), exact-cert checkmarks
Mathtext note: use \\leq (no bare \\le in this matplotlib). No tight_layout
+ aspect('equal') combination (wrecks layout) - manual gridspec instead.
"""
import sympy as sp, numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

rng = np.random.default_rng(5)
w, w1, w2 = sp.symbols("w w1 w2")
d = 5
m = d * (d + 1)
p = sp.expand(-w**d + w**(d - 1) - (3 - sp.Rational(6, m)) * w**2 + (2 - sp.Rational(6, m)) * w)
Ph = sp.expand(sp.integrate(p, w))
pp = sp.diff(p, w)
eq1 = sp.expand(sp.cancel((p.subs(w, w2) - p.subs(w, w1)) / (w2 - w1)))
N2 = sp.expand(Ph.subs(w, w2) - Ph.subs(w, w1) - (w2 - w1) * p.subs(w, w1))
te2 = sp.expand(sp.div(sp.Poly(N2, w1, w2), sp.Poly((w2 - w1)**2, w1, w2))[0].as_expr())

# ---- d=5 exact-ish data ----
flexes = [complex(r) for r in sp.nroots(pp, n=60, maxsteps=200)]
flex_real = sorted([r.real for r in flexes if abs(r.imag) < 1e-40])
print("flexes (all):", [f"{r.real:.10f}{r.imag:+.2e}i" for r in flexes])
print("real flexes:", flex_real)

den5 = 5
E5 = sp.resultant(eq1, (w2 - w1) * te2, w1)
CofQ = sp.Poly(sp.cancel(E5 / pp.subs(w, w2)**2), w2)
cofs = [mp_c for mp_c in CofQ.all_coeffs()]
roots = [complex(r) for r in sp.nroots(CofQ, n=80, maxsteps=300)]
def Pv(z): return complex(sp.N(p.subs(w, z), 50))
def Fv(z): return complex(sp.N(Ph.subs(w, z), 50))
pvs = [Pv(r) for r in roots]
n5 = len(roots)
used = [False] * n5; pairs5 = []
for i in range(n5):
    if used[i]: continue
    best, bj = 1e30, -1
    for j in range(n5):
        if j != i and not used[j] and abs(pvs[i] - pvs[j]) < best:
            best, bj = abs(pvs[i] - pvs[j]), j
    used[i] = used[bj] = True; pairs5.append((roots[i], roots[bj], pvs[i]))

node_lines = []   # (s, r, t1, t2, kind)
for a, b, v in pairs5:
    sv = v                               # p-value at contact (= node s)
    rv = Fv(a) - sv * a                  # node r = Phi(a) - s*a
    kind = ("crun" if abs(a.imag) < 1e-30 and abs(b.imag) < 1e-30
            else "acn" if abs(b - a.conjugate()) < 1e-20 else "cplx")
    if kind == "cplx":
        continue
    s0, r0 = float(sv.real), float(rv.real)
    node_lines.append((s0, r0, a, b, kind))
    print(f"node: s={s0:.12f} r={r0:.12f} t=({a:.6f},{b:.6f}) {kind}")

eq1n = sp.lambdify((w1, w2), eq1, "numpy")
te2n = sp.lambdify((w1, w2), te2, "numpy")

fig = plt.figure(figsize=(14.5, 11.5))
gs = fig.add_gridspec(2, 2, left=0.055, right=0.975, top=0.925, bottom=0.075,
                      wspace=0.21, hspace=0.30)

# ================================================================ (a)
ax = fig.add_subplot(gs[0, 0])
xx = np.linspace(-1.9, 1.45, 700)
yy = np.linspace(-1.9, 1.45, 700)
T1, T2 = np.meshgrid(xx, yy)
Z1 = eq1n(T1, T2); Z2 = te2n(T1, T2)
ax.contour(T1, T2, Z1, levels=[0], colors=["#1a56c4"], linewidths=[1.3])
ax.contour(T1, T2, Z2, levels=[0], colors=["#c41a3a"], linewidths=[1.3])
ax.plot([-1.9, 1.45], [-1.9, 1.45], "k--", lw=1.0, alpha=0.55)
for fr in flex_real:
    ax.plot(fr, fr, "ks", ms=7, mfc="none", mew=1.8, zorder=5)
npts = []
npts_c = []
for a, b, v in pairs5:
    if abs(a.imag) < 1e-30 and abs(b.imag) < 1e-30:
        npts += [(a.real, b.real), (b.real, a.real)]
for (X, Y) in npts:
    ax.plot(X, Y, "o", ms=8, mfc="#ff8c1a", mec="k", mew=1.0, zorder=6)
ax.set_xlim(-1.9, 1.45); ax.set_ylim(-1.9, 1.45)
ax.set_xlabel("$t_1$", fontsize=12); ax.set_ylabel("$t_2$", fontsize=12)
ax.set_title("(a) the incidence chamber at $d=5$:  $\\{eq_1=0\\}\\cap\\{te_2=0\\}$\n"
             "16 = $(d-1)^2$ points over $\\mathbb{C}$: 4 cusps (diag, □) + 12 ordered node-offs (●)",
             fontsize=11)
ax.text(0.03, 0.97, "real points: 2 cusps + 2 crunode offs;\nother 12 intersections complex", transform=ax.transAxes,
        fontsize=9, va="top", ha="left", bbox=dict(fc="white", ec="gray", alpha=0.85))
hs = [Line2D([0], [0], color="#1a56c4", lw=1.5, label="$eq_1$: $p(t_1)=p(t_2)$  (slope $-1$ at flex)"),
      Line2D([0], [0], color="#c41a3a", lw=1.5, label="$te_2$: second divided diff (slope $-2$)"),
      Line2D([0], [0], color="k", ls="--", lw=1.2, label="diagonal $t_1=t_2$ (slope $+1$)")]
ax.legend(handles=hs, fontsize=8.3, loc="lower left", framealpha=0.9)
ax.set_aspect("equal")

# ================================================================ (b)
tau = flex_real[0]
print("zoom flex tau =", tau)
axb = fig.add_subplot(gs[0, 1])
eps = 0.006
xb = np.linspace(tau - eps, tau + eps, 400)
yb = np.linspace(tau - eps, tau + eps, 400)
TB1, TB2 = np.meshgrid(xb, yb)
axb.contour(TB1, TB2, eq1n(TB1, TB2), levels=[0], colors=["#1a56c4"], linewidths=[1.6])
axb.contour(TB1, TB2, te2n(TB1, TB2), levels=[0], colors=["#c41a3a"], linewidths=[1.6])
axb.plot([tau - eps, tau + eps], [tau - eps, tau + eps], "k--", lw=1.2)
xxs = np.linspace(tau - eps, tau + eps, 3)
axb.plot(xxs, tau - (xxs - tau), color="#1a56c4", lw=0.8, ls=":", alpha=0.9)
axb.plot(xxs, tau - 2 * (xxs - tau), color="#c41a3a", lw=0.8, ls=":", alpha=0.9)
axb.plot(tau, tau, "ks", ms=9, mfc="none", mew=2.2, zorder=6)
axb.text(tau - 0.90 * eps, tau + 0.84 * eps, "slope $-1$\n$eq_1$: $3c_3(u+x)+O(2)$", fontsize=9.5,
         color="#1a56c4", rotation=-45, ha="center", va="center")
axb.text(tau + 0.40 * eps, tau - 0.84 * eps, "slope $-2$\n$te_2$: $c_3(u+2x)+O(2)$", fontsize=9.5,
         color="#c41a3a", rotation=-63, ha="center", va="center")
axb.text(tau + 0.50 * eps, tau + 0.68 * eps, "diag $+1$\n$(u-x)$", fontsize=9.5, color="k",
         rotation=45, ha="center", va="center")
axb.set_xlim(tau - eps, tau + eps); axb.set_ylim(tau - eps, tau + eps)
axb.set_xlabel("$t_1$", fontsize=12); axb.set_ylabel("$t_2$", fontsize=12)
c3 = float(sp.N(sp.diff(p, w, 2).subs(w, sp.Float(tau, 40)) / 6, 12))
axb.set_title(f"(b) universal germ at flex $\\tau={tau:.6f}$:  $c_3=p''(\\tau)/6={c3:.4f}\\neq0$\n"
              "three branches transverse $\\Rightarrow$ the cusp toll is exactly 2 (all $d$, proven)",
              fontsize=11)
axb.set_aspect("equal")

# ================================================================ (c)
axc = fig.add_subplot(gs[1, 0])
wg = np.linspace(-1.95, 1.22, 1500)
Phf = sp.lambdify(w, Ph, "numpy")
axc.plot(wg, Phf(wg), color="#0a6e2e", lw=1.7, label="$\\Gamma_5: \\ r=\\Phi_5(w)$")
for fr in flex_real:
    s0 = float(sp.N(p.subs(w, sp.Float(fr, 40)), 25))
    r0 = float(sp.N(Ph.subs(w, sp.Float(fr, 40)), 25))
    xs = np.linspace(-1.85, 1.22, 4)
    axc.plot(xs, r0 + s0 * (xs - fr), color="#c41a3a", lw=1.2, alpha=0.85)
    axc.plot(fr, r0, "s", mfc="none", mec="#c41a3a", ms=8, mew=1.8)
lab_map = {"crun": ("#ff8c1a", "-", 2.2), "acn": ("#7a1ac4", "--", 1.4)}
for s0, r0, a, b, kind in node_lines:
    xs = np.linspace(-1.85, 1.22, 4)
    col, lsc, lwc = lab_map[kind]
    axc.plot(xs, r0 + s0 * xs, color=col, ls=lsc, lw=lwc, alpha=0.95)
    if kind == "crun":
        for tt in (a.real, b.real):
            axc.plot(tt, r0 + s0 * tt, "o", ms=7, mfc=col, mec="k")
axc.set_xlim(-1.95, 1.22); axc.set_ylim(-2.6, 2.2)
axc.axhline(0, color="gray", lw=0.6, alpha=0.4)
axc.set_xlabel("$w$", fontsize=12); axc.set_ylabel("$r$", fontsize=12)
axc.set_title("(c) self-portrait of the graph of $\\Phi_5$: 2 real flex tangents (red),\n"
              "crunode bitangent (orange, 2 real contacts), acnode line (purple, complex contacts)",
              fontsize=11)
hc = [Line2D([0], [0], color="#0a6e2e", lw=1.7, label="$\\Phi_5(w) = -w^6/6+w^5/5-\\frac{14}{15}w^3+\\frac{9}{10}w^2$"),
      Line2D([0], [0], color="#c41a3a", lw=1.3, label="flex tangent -> wall cusp"),
      Line2D([0], [0], color="#ff8c1a", lw=2.2, label="crunode bitangent -> wall node"),
      Line2D([0], [0], color="#7a1ac4", lw=1.4, ls="--", label="acnode bitangent (real line, $\\mathbb{C}$ contacts)")]
axc.legend(handles=hc, fontsize=8.3, loc="lower left", framealpha=0.92)

# ================================================================ (d)
axd = fig.add_subplot(gs[1, 1])
col = json.load(open("budgetlaw_stage2.json"))["column"]
ds = [dd for dd in range(4, 13)]
cuspT = [2 * (dd - 1) for dd in ds]
nodeT = [(dd - 1) * (dd - 2) for dd in ds]
Ks = [col[str(dd)]["K"] for dd in ds]
xpos = np.arange(len(ds))
axd.bar(xpos, cuspT, 0.62, color="#ff8c1a", label="cusp toll  $2(d-1)$")
axd.bar(xpos, nodeT, 0.62, bottom=cuspT, color="#1a56c4", alpha=0.85, label="node toll  $(d-1)(d-2)$  (= 2 orientations $\\times$ #(pairs))")
for i, dd in enumerate(ds):
    tot = cuspT[i] + nodeT[i]
    axd.text(xpos[i], tot + 2.0, f"{tot}", ha="center", fontsize=9, fontweight="bold")
    axd.text(xpos[i], -9.5, f"den={col[str(dd)]['den']}", ha="center", fontsize=7.5, color="dimgray")
    ok = col[str(dd)]["sqfree"] and col[str(dd)]["coprime"]
    axd.text(xpos[i], tot + 8.5, "✓✓" if ok else "✗", ha="center", fontsize=10,
             color="#0a6e2e" if ok else "red")
axd.set_xticks(xpos)
axd.set_xticklabels([f"d={dd}\nn={dd+1}" for dd in ds], fontsize=8.5)
axd.set_ylabel("deg $E_d(t_2)$  (budget $d(d-1)$)", fontsize=11)
axd.set_ylim(-14, 155)
axd2 = axd.twinx()
axd2.plot(xpos, Ks, "D", color="#c41a3a", ms=7, mec="k", label="$K=\\mathrm{den}(p'_d)^2$")
axd2.set_yscale("log")
axd2.set_ylabel("$K$ (log scale)", fontsize=11, color="#c41a3a")
axd2.tick_params(axis="y", colors="#c41a3a")
for i, K in enumerate(Ks):
    axd2.annotate(f"{K}", xy=(xpos[i], K), xytext=(xpos[i] + 0.13, K * 1.4), fontsize=7.5, color="#c41a3a")
axd.set_title("(d) evidence column: $E_d=p_d'^2\\,\\mathrm{Cof}_d$ exact, $d=4..13$;\n"
              "✓✓ = $\\mathrm{Cof}$ squarefree ∧ coprime to $p'_d$ (T1 $\\wedge$ T2), exact gcd, $d=4..12$",
              fontsize=11)
hd, ld = axd.get_legend_handles_labels()
hd2, ld2 = axd2.get_legend_handles_labels()
axd.legend(hd + hd2, ld + ld2, fontsize=8.6, loc="upper left", framealpha=0.92)

fig.suptitle("NOTE 21 — THE BUDGET LAW:  $E_d(t_2)=\\mathrm{Res}_{t_1}(eq_1,eq_2)=p_d'(t_2)^2\\cdot\\mathrm{Cof}_d(t_2)$   "
             "·   every wall singularity pays multiplicity exactly 2", fontsize=13, y=0.985)
fig.savefig("budgetlaw_figure.png", dpi=150)
print("figure saved")
