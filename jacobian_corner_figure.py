"""NOTE 15 figure: the corner's address, four panels."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sympy as sp
from sympy import symbols, expand, integrate, Rational as R

w = symbols("w")
def seed(d):
    return expand(2*w - 3*w**2 + w*(1-w)*(w**sp.Integer(d-2) - R(6, d*(d+1))))

fig, axes = plt.subplots(2, 2, figsize=(14.6, 12.6))
plt.subplots_adjust(hspace=0.30, wspace=0.24, left=0.065, right=0.975, top=0.93, bottom=0.06)

# ---------------- panel A: the d=9 wall w/ full entourage ----------------
ax = axes[0, 0]
d = 9
p = seed(d); Phi = integrate(p, w)
pf = sp.lambdify(w, p, "numpy"); phif = sp.lambdify(w, Phi, "numpy")
T = np.concatenate([np.linspace(-0.998, 2.6, 2400), np.linspace(-2.6, -1.002, 2400)])
Sv = pf(T); Rv = T*pf(T) - phif(T)
# break the polyline at huge jumps (poles of the sampling)
pts = np.array([Sv, Rv])
ax.plot(Sv, Rv, '-', color='0.35', lw=0.9, zorder=2, label="wall $t\\mapsto(p(t),\\tau(t))$, d=9")
ax.plot([-4.3, 3.4], [-4.3, 3.4], ':', color='firebrick', lw=1.1, zorder=1)
ax.text(1.55, 1.2, "diagonal $r=s$", color='firebrick', fontsize=9, rotation=32)
j2 = json.load(open("jcorner_basins.json"))
def scr(x): return complex(str(x).replace(" ", ""))
for nd in j2["basins"]["9"]:
    scv, rcv = scr(nd[1]).real, scr(nd[2]).real
    if nd[0] == "CRUNODE":
        ax.plot(scv, rcv, 's', ms=11, mfc='none', mec='navy', mew=2.2, zorder=6)
        ax.annotate(f"CRUNODE (corner)\n({scv:.4f}, {rcv:.4f})", (scv, rcv),
                    textcoords="offset points", xytext=(-118, 26), fontsize=9, color='navy')
    else:
        ax.plot(scv, rcv, 'o', ms=7, mfc='none', mec='darkorange', mew=1.6, zorder=5)
# markers: pin, shadow, cusps, t=-1, origin
import mpmath as mp
s1 = json.load(open('jcorner_stage1.json'))
reload = s1
tv9 = [mp.mpf(q['t']) for q in reload['diag']['9'] if mp.mpf(q['t']) < -1.0001][0]
s9 = mp.mpf([q['s'] for q in reload['diag']['9'] if mp.mpf(q['t']) < -1.0001][0])
ax.plot(-1, -1, 'P', ms=13, color='green', zorder=7)
ax.annotate("pin $t=1$: $(-1,-1)$\non EVERY tower wall", (-1, -1), textcoords="offset points",
            xytext=(26, -34), fontsize=9, color='green')
ax.plot(float(s9), float(s9), 'X', ms=11, color='purple', zorder=7)
ax.annotate(f"diagonal shadow\n$s^*={float(s9):.5f}$\n$t^*={float(tv9):.5f}$", (float(s9), float(s9)),
            textcoords="offset points", xytext=(-132, -48), fontsize=9, color='purple')
# cusps
for tgv, lab, col, offs in [(1/3, "ghost cusp image\n$\\to(1/3,1/27)$", 'crimson', (-130, 6)),
                               (float(mp.mpf(s1['races']['9']['t_left'])), "left cusp image\n$\\to(-5,3)$", 'crimson', (-60, 42))]:
    s_c = float(pf(tgv)); r_c = float(tgv*pf(tgv) - phif(tgv))
    ax.plot(s_c, r_c, '*', ms=13, color=col, zorder=6)
    ax.annotate(lab, (s_c, r_c), textcoords="offset points", xytext=offs, fontsize=8.6, color=col)
# t=-1 image
s_m1 = float(pf(-1)); r_m1 = float(-1*pf(-1) - phif(-1))
ax.plot(s_m1, r_m1, 'd', ms=10, mfc='none', mec='teal', mew=2, zorder=6)
ax.annotate("$t=-1$ image", (s_m1, r_m1), textcoords="offset points", xytext=(-4, 22), fontsize=8.6, color='teal')
ax.plot(0, 0, 'h', ms=10, mfc='none', mec='dimgray', mew=2, zorder=6)
ax.annotate("hole node (0,0)", (0, 0), textcoords="offset points", xytext=(12, 12), fontsize=8.6, color='dimgray')
ax.set_xlim(-5.4, 3.1); ax.set_ylim(-2.6, 3.6)
ax.set_xlabel("s"); ax.set_ylabel("r")
ax.set_title("(A) the wall of chamber n=10 (d=9) and its entourage", fontsize=11.5)
ax.grid(alpha=0.25)

# ---------------- panel B: corner zoom ----------------
ax = axes[0, 1]
T2 = np.linspace(-2.2, 2.2, 5000)
Sv = pf(T2); Rv = T2*pf(T2) - phif(T2)
ax.plot(Sv, Rv, '-', color='0.35', lw=1.1, label="wall d=9")
ax.plot([-1.18, -0.70], [-1.18, -0.70], ':', color='firebrick', lw=1.3)
ax.plot(-1, -1, 'P', ms=14, color='green', zorder=7)
ax.annotate("pin: diagonal is TANGENT\nto the wall here (slope $t=1$)", (-1, -1),
            textcoords="offset points", xytext=(-160, 44), fontsize=9, color='green')
ax.plot(float(s9), float(s9), 'X', ms=12, color='purple', zorder=7)
ax.annotate("shadow", (float(s9), float(s9)), textcoords="offset points", xytext=(-64, -4), fontsize=9, color='purple')
crun = [nd for nd in j2["basins"]["9"] if nd[0] == "CRUNODE"][0]
cs, cr_ = scr(crun[1]).real, scr(crun[2]).real
ax.plot(cs, cr_, 's', ms=10, mfc='none', mec='navy', mew=2, zorder=7)
ax.annotate("crunode", (cs, cr_), textcoords="offset points", xytext=(8, 12), fontsize=9, color='navy')
t1c, t2c = scr(crun[3]).real, scr(crun[4]).real
ax.annotate(f"contacts $t_1={t1c:+.5f},\ t_2={t2c:+.5f}$\nboth map HERE (node!)\n$|t_2-t^*|$ ~ $0.0066/d^3$",
            (cs, cr_), textcoords="offset points", xytext=(-165, -64), fontsize=8.4, color='darkorange')
ax.set_xlim(-1.155, -0.80); ax.set_ylim(-1.155, -0.80)
ax.set_xlabel("s"); ax.set_ylabel("r")
ax.set_title("(B) corner zoom: pin $=$ tangency $+$ shadow $+$ crunode coalesce", fontsize=11.5)
ax.grid(alpha=0.25)

# ---------------- panel C: the races ----------------
ax = axes[1, 0]
st1 = json.load(open('jcorner_stage1.json'))
Ds_sh, U_sh, G5, Gc, Lc = [], [], [], [], []
Dn_ = sorted(int(k) for k in st1['races'].keys())
for d_ in Dn_:
    if d_ < 5: continue
    r = st1['races'][str(d_)]
    u_left = mp.mpf(r['u'])
    Lc.append(u_left * (d_-2) / mp.log((d_-2)/4))
# exact shadow set from stage json
fin = json.load(open('jcorner_final.json'))
# recompute shadow locally (cheap for odd d<=45)
import numpy as np
for d_ in sorted([3] + list(range(5, 46, 2))):
    pp = seed(d_); Phi_ = integrate(pp, w)
    G = expand(Phi_ - (w-1)*pp)
    cfl = [float(c) for c in sp.Poly(sp.expand(G*(d_*(d_+1))), w).all_coeffs()]
    root = None
    for z in np.roots(cfl):
        if abs(z.imag) < 1e-7 and z.real < -1.0001 and (root is None or z.real > root):
            root = z.real
    if d_ == 3:
        continue
    Ds_sh.append(d_)
    U_sh.append(abs(root+1)*(d_-2))
    pass
Ln2 = float(np.log(2))
ax.plot(Ds_sh, U_sh, 'o-', ms=5, lw=1.2, color='navy', label="shadow: $U(d)\\,(d-2)\\to\\ln 2$")
ax.axhline(Ln2, color='navy', ls=':', lw=1)
ax.text(45.5, Ln2, "$\\ln 2$", color='navy', va='center', fontsize=10)
dsr = [d_ for d_ in Dn_ if d_ >= 5]
R_ = [float(mp.mpf(st1['races'][str(d_)]['u'])*(d_-2)/mp.log((d_-2)/4)) for d_ in dsr]
ax.plot(dsr, R_, 's-', ms=5, lw=1.2, color='darkred', label="left cusp: $u(d)(d-2)/\\ln\\frac{d-2}{4}\\to 1$")
ax.axhline(1, color='darkred', ls=':', lw=1)
dg_ = [d_ for d_ in range(5, 42, 2)]
gser = [float(mp.mpf(st1['races'][str(d_)]['t_ghost'])) for d_ in dg_]
grc = [float((mp.mpf(1)/3 - t_)* (3*((d_)*(d_+1)-2))) for d_, t_ in zip(dg_, gser)]
ax.plot(dg_, grc, 'd-', ms=5, lw=1.2, color='darkgreen', label="ghost: $(1/3-t_g)\\,3[d(d+1)-2]\\to 1$")
sgd = [d_ for d_ in sorted([5] + list(range(5, 46, 2)))]
sigv = []
for d_ in sgd:
    pp = seed(d_); Phi_ = integrate(pp, w)
    G = expand(Phi_ - (w-1)*pp)
    cfl = [float(c) for c in sp.Poly(sp.expand(G*(d_*(d_+1))), w).all_coeffs()]
    root = None
    for z in np.roots(cfl):
        if abs(z.imag) < 1e-7 and z.real < -1.0001 and (root is None or z.real > root):
            root = z.real
    sv = float(pp.subs(w, sp.Float(root, 25)))
    sigv.append((sv+1)*(d_-2))
ax.plot(sgd, sigv, '^-', ms=5, lw=1.2, color='purple', label="shadow s-cord: $(s^*+1)(d-2)\\to 2-2\\ln 2$")
ax.axhline(2-2*Ln2, color='purple', ls=':', lw=1)
ax.text(45.5, 2-2*Ln2, "$2-2\\ln 2$", color='purple', va='center', fontsize=10)
ax.set_xlabel("d (odd)"); ax.set_ylabel("race profile")
ax.set_title("(C) the four races of the corner basin", fontsize=11.5)
ax.legend(fontsize=8.6, loc='upper right')
ax.grid(alpha=0.25); ax.set_ylim(0, 2.3)

# ---------------- panel D: d=5 exact count-map teardrop ----------------
ax = axes[1, 1]
z = np.load("jcorner_map_d5.npz")
S_, R_, G_ = z["S"], z["R"], z["G"]
cmap = matplotlib.colors.ListedColormap(['#ffd54f', '#9fb7d9', '#bbdefb', '#1e3a5f', '#90caf9'])
imv = np.clip(G_, 0, 4)
ax.imshow(imv, origin='lower', extent=[S_[0], S_[-1], R_[0], R_[-1]], cmap=cmap, aspect='auto')
ax.plot([-1.15, -0.35], [-1.15, -0.35], ':', color='firebrick', lw=1.1)
j5 = [nd for nd in j2["basins"]["5"] if nd[0] == "CRUNODE"][0]
ax.plot(scr(j5[1]).real, scr(j5[2]).real, 's', ms=12, mfc='none', mec='navy', mew=2.2)
ax.annotate(f"crunode ({scr(j5[1]).real:.4f},{scr(j5[2]).real:.4f})\n= tip of the missed teardrop",
            (scr(j5[1]).real, scr(j5[2]).real), textcoords="offset points", xytext=(-150, -38),
            fontsize=9, color='navy')
s5 = [q for q in reload['diag']['5'] if mp.mpf(q['t']) < -1.0001][0]
ax.plot(float(mp.mpf(s5['s'])), float(mp.mpf(s5['s'])), 'X', ms=11, color='purple')
ax.annotate("shadow (on diagonal)", (float(mp.mpf(s5['s'])), float(mp.mpf(s5['s']))),
            textcoords="offset points", xytext=(10, 10), fontsize=9, color='purple')
ax.plot(-1, -1, 'P', ms=12, color='green')
from matplotlib.lines import Line2D
leg = [Line2D([0],[0], marker='s', color='w', mfc='#ffd54f', ms=11, label='0 real fibers (MISSED)'),
       Line2D([0],[0], marker='s', color='w', mfc='#1e3a5f', ms=11, label='4 real fibers'),
       Line2D([0],[0], marker='s', color='w', mfc='#bbdefb', ms=11, label='2 real fibers')]
ax.legend(handles=leg, fontsize=8.6, loc='upper right')
ax.set_xlabel("s"); ax.set_ylabel("r")
ax.set_title("(D) d=5: the cone is a teardrop hung from its crunode tip\n(exact count-map, $1e^{{\\dots}}$ audited machine Sturm)", fontsize=11.5)

fig.suptitle("THE CORNER'S ADDRESS — the miss$(-1,-1)$ basin's algebraic census and its migration races", fontsize=14, y=0.975)
fig.savefig("jcorner_figure.png", dpi=150)
print("saved jcorner_figure.png")
