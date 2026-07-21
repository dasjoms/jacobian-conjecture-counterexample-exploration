"""Figure for note 6: the escape atlas of F4."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C5, C4, C3, C2 = -1/5, 1/4, -9/10, 17/20
def p4v(w): return -w**4 + w**3 - 2.7*w**2 + 1.7*w
def Phiv(w): return -w**5/5 + w**4/4 - 0.9*w**3 + 0.85*w**2
D4_terms = [(-20000000,4,0),(20000000,3,1),(18900000,3,0),(-57100000,2,2),(75376000,2,1),
            (-46613461,2,0),(38210000,1,3),(-67741050,1,2),(45893931,1,1),(-6657115,1,0),
            (-8192000,0,5),(17058675,0,4),(-12278715,0,3),(1957975,0,2)]
def discv(s, r):
    return sum(c * r**rr * s**ss for c, rr, ss in D4_terms)

def track_trajs(path, steps=2200):
    ts = np.linspace(0, 2*np.pi, steps, endpoint=False)
    pts = [path(t) for t in ts]
    def greedy(prev, nxt):
        used=[False]*5; out=[None]*5
        for a in range(5):
            bd,bj=1e30,0
            for b in range(5):
                if not used[b] and abs(prev[a]-nxt[b])<bd: bd,bj=abs(prev[a]-nxt[b]),b
            used[bj]=True; out[a]=nxt[bj]
        return out
    cur = list(np.roots([C5,C4,C3,C2,-pts[0][0],pts[0][1]]))
    tr = [[c] for c in cur]
    for i in range(1, len(pts)):
        cur = greedy(cur, list(np.roots([C5,C4,C3,C2,-pts[i][0],pts[i][1]])))
        for k in range(5): tr[k].append(cur[k])
    return tr

fig, axes = plt.subplots(2, 2, figsize=(13, 11))

# ---- (a) atlas ----
ax = axes[0][0]
tt = np.linspace(-2.4, 3.4, 3000)
ax.plot(p4v(tt), tt*p4v(tt)-Phiv(tt), color="#1f77b4", lw=1.8,
        label="real tangent lines $t\\mapsto(p(t),\\, t p(t)-\\Phi(t))$")
gs = np.linspace(-2.6, 2.6, 500); gr = np.linspace(-1.4, 2.1, 400)
S, Rg = np.meshgrid(gs, gr)
DD = discv(S, Rg)
ax.contour(S, Rg, DD, levels=[0], colors="#bbbbbb", linewidths=0.8)
ax.contour(S, Rg, DD, levels=[-1e5,1e5], colors="none")
marks = [
    (0, 0, "o", "k", "$(0,0)$: all of $\{C=0\}$", (8,-10)),
    (-1, -1, "s", "seagreen", "fold $(-1,-1) = t{=}1$", (10,-12)),
    (0.2921225227866996, 0.0340042595568076877, "X", "crimson", "real cusp $(3)$", (12,-4)),
    (0.9841655956599125, 0.67623935, "v", "darkviolet", "acnode $(2,2)$, complex contacts", (-30,12)),
]
for sx, ry, mk, cl, lb, off in marks:
    ax.plot([sx],[ry],mk,color=cl,ms=9 if mk!="X" else 11)
    ax.annotate(lb, xy=(sx,ry), textcoords="offset points", xytext=off, fontsize=9, color=cl)
for (cx, cy, rr) in [(-1,-1,0.35),(0.29212,0.034,0.28),(0.98417,0.67624,0.3)]:
    th = np.linspace(0,2*np.pi,100)
    ax.plot(cx+rr*np.cos(th), cy+rr*np.sin(th), "g--", lw=1, alpha=0.7)
ax.annotate("monodromy loops\n(schematic radii)", xy=(-1+0.35,-0.75), fontsize=8, color="g")
ax.set_xlim(-2.6,2.6); ax.set_ylim(-1.5,2.2)
ax.set_xlabel("s = BC"); ax.set_ylabel("r = AC$^2$")
ax.set_title("(a) THE WALL of $F_4$: $D_4(s,r)=0$ (grey full real trace; blue = tangent lines)\n"
             "rational quintic, 3 cusps + 3 bitangent nodes = maximally singular", fontsize=10.5)
ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="lower left")

# ---- (b) fold loop trajectories ----
ax = axes[0][1]
tr = track_trajs(lambda t: (-1 + 0.03*np.exp(1j*t)/np.sqrt(2), -1 + 1j*0.03*np.exp(1j*t)/np.sqrt(2)))
cols = plt.cm.tab10.colors
for k, traj in enumerate(tr):
    A = np.array(traj)
    ax.plot(A.real, A.imag, color=cols[k], lw=1.2)
    ax.plot([A[0].real],[A[0].imag],"o",color=cols[k],ms=7)
    ax.plot([A[-1].real],[A[-1].imag],"^",color=cols[k],ms=7)
ax.set_xlabel("Re w"); ax.set_ylabel("Im w")
ax.set_title("(b) fold loop around $(-1,-1)$: one root-pair exchanges\n(tiny ellipse, right); the other 3 hardly move", fontsize=10.5)
ax.grid(alpha=0.25)

# ---- (c) cusp loop trajectories ----
ax = axes[1][0]
tr = track_trajs(lambda t: (0.2921225227866996 + 0.008*np.exp(1j*t)/np.sqrt(2),
                            0.0340042595568076877 + 1j*0.008*np.exp(1j*t)/np.sqrt(2)))
for k, traj in enumerate(tr):
    A = np.array(traj)
    ax.plot(A.real, A.imag, color=cols[k], lw=1.2)
    ax.plot([A[0].real],[A[0].imag],"o",color=cols[k],ms=7)
    ax.plot([A[-1].real],[A[-1].imag],"^",color=cols[k],ms=7)
ax.set_xlabel("Re w"); ax.set_ylabel("Im w")
ax.set_title("(c) cusp loop: 3 roots braid into one long 3-cycle\n(central ellipse); two at the corners stay put", fontsize=10.5)
ax.grid(alpha=0.25)

# ---- (d) escape rate at fold ----
ax = axes[1][1]
deltas = np.array([1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8])
mins = []
for delta in deltas:
    s, r = -1.0, -1.0 + delta
    rts = np.roots([C5,C4,C3,C2,-s,r])
    gam = np.abs(s - p4v(rts))
    gam.sort()
    mins.append(gam[0])
mins = np.array(mins)
ax.loglog(deltas, mins, "o-", color="#1f77b4", label="$\\min|\\gamma|$ over 5 roots")
ax.loglog(deltas, 3.07*np.sqrt(deltas), "--", color="crimson", label="$3.07\\,\\delta^{1/2}$")
ax.set_xlabel("wall distance $\\delta$"); ax.set_ylabel("$|\\gamma|$ of escaping sheet")
ax.set_title("(d) fold escape rate: $\\gamma \\sim 3.07\\sqrt{\\delta}$, so\n$|x| = |C/\\gamma| \\sim \\delta^{-1/2}$: the sheets blow up with slope $-1/2$", fontsize=10.5)
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig("/home/user/escape_atlas.png", dpi=140)
print("figure saved")
