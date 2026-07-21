"""
Monodromy of the 3-sheeted covering: track the roots of the fiber cubic
w^3 - w^2 + P w - Q = 0 as (P,Q) travels around loops in parameter space.
Expectation: loops not linked with the wall curve Gamma give the identity;
loops crossing one branch give a transposition; loops around the CUSP may
give a 3-cycle. (Together: monodromy = S_3.)
"""
import numpy as np

def fiber_roots(P, Q):
    return np.roots([1, -1, P, -Q])          # w^3 - w^2 + Pw - Q

def track(c, r, n=1500):
    """loop: (P,Q) = c + r*(cos t, sin t); returns end permutation + trajectories"""
    ts = np.linspace(0, 2 * np.pi, n)
    prev = fiber_roots(c[0] + r * np.cos(0), c[1] + r * np.sin(0))
    start = prev.copy()
    traj = [prev]
    for t in ts[1:]:
        cur = fiber_roots(c[0] + r * np.cos(t), c[1] + r * np.sin(t))
        # greedy matching to previous
        used = set(); new = np.zeros(3, complex)
        for j, pr in enumerate(prev):
            k = min((i for i in range(3) if i not in used), key=lambda i: abs(cur[i] - pr))
            used.add(k); new[j] = cur[k]
        prev = new
        traj.append(prev)
    perm = [int(np.argmin([abs(start[j] - prev[i]) for j in range(3)])) for i in range(3)]
    return perm, np.array(traj)

def R(P, Q):   # wall curve
    return 4 * P**3 - P**2 - 18 * P * Q + 27 * Q**2 + 4 * Q

tests = {
    "A: big loop, center (0.1,-0.05), r=1.4": ((0.1, -0.05), 1.4),
    "B: center (0.1,-0.05), r=0.5":          ((0.1, -0.05), 0.5),
    "C: small loop around cusp (1/3,1/27), r=0.08": ((1/3, 1/27), 0.08),
    "D: tiny loop at generic point (1.5, 0.2), r=0.05": ((1.5, 0.2), 0.05),
}
for name, (c, r) in tests.items():
    perm, traj = track(np.array(c, complex) if isinstance(c, tuple) else c, r)
    print(f"{name}: end permutation = {perm}  "
          f"({'identity' if perm == [0,1,2] else 'nontrivial: ' + str(perm)})")
    # crossing diagnostics: sign changes of R along the loop
    import numpy as _np
    ts = _np.linspace(0, 2*_np.pi, 4000, endpoint=False)
    vals = _np.array([R(c[0] + r*_np.cos(t), c[1] + r*_np.sin(t)) for t in ts])
    crossings = int(_np.sum(vals[:-1] * vals[1:] < 0))
    print(f"    wall-curve crossings along loop: {crossings}")

# save trajectory figure for loops B and C
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, key in zip(axes, ["B: center (0.1,-0.05), r=0.5", "C: small loop around cusp (1/3,1/27), r=0.08"]):
    c, r = tests[key]; _, traj = track(c, r)
    for j in range(3):
        ax.plot(traj[:, j].real, traj[:, j].imag, ".", ms=2, label=f"sheet {j+1}")
    midd = len(traj) // 2
    ax.plot(traj[0, :].real, traj[0, :].imag, "ko", ms=5, label="_start")
    ax.plot(traj[-1, :].real, traj[-1, :].imag, "k^", ms=6, label="_end")
    for j in range(3):
        ax.annotate(str(j+1), (traj[0, j].real, traj[0, j].imag), fontsize=8)
    ax.set_title(key.split(":")[0] + ": root trajectories in w-plane")
    ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig("/home/user/monodromy.png", dpi=130)
print("figure: monodromy.png")
