import sympy as sp
from sympy import symbols, diff, integrate, expand, Rational as R
w = symbols("w")
p4 = -w**4 + w**3 - R(27,10)*w**2 + R(17,10)*w
Phi4 = expand(integrate(p4, w))
sext = 6400*w**6 - 9600*w**5 + 33280*w**4 - 34600*w**3 + 47269*w**2 - 29679*w + 13550
rts = [complex(v) for v in sp.nroots(sp.Poly(sext, w), n=40, maxsteps=600)]
print("6 contact points:", [f"{v:.6f}" for v in rts])
info = []
for v in rts:
    sv = complex(p4.subs(w, v)); rv = complex(v*p4.subs(w, v) - Phi4.subs(w, v))
    info.append((v, sv, rv))
# pair by matching (s, r)
used = [False]*6
for i in range(6):
    if used[i]: continue
    for j in range(i+1, 6):
        if used[j]: continue
        vi, si, ri = info[i]; vj, sj, rj = info[j]
        if abs(si-sj) < 1e-6 and abs(ri-rj) < 1e-6:
            used[i] = used[j] = True
            real_line = abs(si.imag) < 1e-9 and abs(ri.imag) < 1e-9
            real_contacts = abs(vi.imag) < 1e-9 and abs(vj.imag) < 1e-9
            print(f"  BITANGENT: w1={vi:.6f} w2={vj:.6f}  (s,r)=({si:.8f},{ri:.8f})"
                  f"  real-line={real_line} real-contacts={real_contacts}")
            break
for i in range(6):
    if not used[i]:
        vi, si, ri = info[i]
        print(f"  UNPAIRED: w={vi:.6f} (s,r)=({si:.8f},{ri:.8f})")
