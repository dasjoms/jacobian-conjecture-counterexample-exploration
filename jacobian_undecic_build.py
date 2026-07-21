"""Assemble NOTE 14 (The Whisker Chamber) with the figure embedded."""
import json, base64

fig = base64.b64encode(open("atlas10_figure.png", "rb").read()).decode()
A  = json.load(open("atlas10_stageA.json"))
B  = json.load(open("atlas10_bitangents.json"))
cen = json.load(open("atlas10_realcensus.json"))

bd = {n[4]: n for n in B["nodes"] if n[4] != "complex"}
acr = [(float(n[2].split()[0][1:]), float(n[3].split()[0][1:])) for n in bd.values()]

md = f"""# JACOBIAN LAB — Note 14: **The Whisker Chamber** 🚂🌬️
### Fiber eleven: 64 terms, K = 3025, zero crunodes, four acnodes, and $S_{{11}}$ — every tower law examined on one chamber, simultaneously

*Date: 2026-07-21. Instruments: sympy 1.14 (resultants, Gröbner, gcds, Sturm-exact), mpmath 110–120
digits, numpy (200k census, monodromy tracking). Locks frozen in `jacobian_atlas10_{{1,2,properness,
monodromy,theorem}}.py` before each computation; data in `atlas10_*.json`; figure `atlas10_figure.png`.*

---

## 0. The recipe at d = 10

Explainer seed: $p_d(w) = 2w-3w^2 + w(1-w)\\big(w^{{d-2}} - \\tfrac{{6}}{{d(d+1)}}\\big)$, here

$$p_{{10}}(w) = -w^{{10}} + w^9 - \\tfrac{{162}}{{55}}w^2 + \\tfrac{{107}}{{55}}w,\\qquad
\\Phi_{{11}}(w) = -\\tfrac{{w^{{11}}}}{{11}} + \\tfrac{{w^{{10}}}}{{10}} - \\tfrac{{54}}{{55}}w^3 + \\tfrac{{107}}{{110}}w^2 .$$

Pins $p(1) = -1$, $\\Phi(1) = 0$ (checked exactly), **$\\kappa = p'(1) = -272/55$**, recipe
**$a = -(1+\\kappa)/(2+\\kappa) = -217/162$**, and $\\det JF = bc = 1$ at **5/5** exact rational points
(d = 10, closing the det queue: chambers d = 6, 7, 9, 10 all pointwise-exact now).

Odd fiber (n = 11) ⟹ the missed set is a **whisker** (measure zero), not a cone — the parity-mechanism
of note 11 made flesh one chamber later.

---

## 1. The wall D10 — term law becomes exact

Resultant of $h = 110(\\Phi_{{11}} - sw + r)$ against $h_w$:

* **degree 11**, **irreducible** over ℚ, parametrized by $t \\mapsto (p_{{10}}(t),\\ \\tau(t) = tp_{{10}}(t)-\\Phi_{{11}}(t))$ ✓
* **EXACTLY 64 TERMS** — note 12's lock $66 - 2 = 64$ (term law $n(n+1)/2 - 2$ at $n=11$) ✓
* $D_{{10}}(0,0) = 0$ and $s^2 \\mid D_{{10}}(s,0)$ (holes 1+2, from $w^2 \\mid \\Phi_{{11}}$) ✓
* magnitude law: $|[s^{{11}}]\\,\\mathrm{{res}}_{{raw}}| = 110^{{21}}\\cdot 10^{{10}}\\cdot(1/11)^{{10}}$ — ratio **1** ✓ ($L = 110$)

And tonight the law sharpens into its **exact form**, verified symbolically on both stored walls:

$$\\operatorname{{supp}}(D_n) \\;=\\; \\mathrm{{cone}}(n) \\setminus \\big\\{{\\textstyle{{(0,0),\\, (1,0),\\, (0,n-2),\\, (0,n)}}\\big\\}} \\;\\cup\\; \\big\\{{(n,0)\\big\\}},\\qquad
\\#\\operatorname{{supp}} = \\frac{{n(n+1)}}2 + 1 - 3$$

where $\\mathrm{{cone}} = \\{{(n-1)j + ni \\le n(n-1)\\}}$ ($n(n+1)/2 + 1$ lattice points). The fourth "hole"
$r^{{n}}$ is **fictitious** — the resultant's $r$-degree is $n-1$ by construction, so that lattice point
can never be populated — and it is precisely balanced by $s^n$ reaching **one weight step outside the cone**
(the same $s^n$ whose coefficient the magnitude law governs). The three genuine holes are const, $s$, and
$r^{{n-2}}$ — the third hole persisting as the $\\beta$-theorem predicted.

* fingerprints: $[s^{{11}}]D_{{10}} = -2^{{20}}\\,5^{{28}}\\,11^{{9}}$ — sign $(-1)^n$ holds 9/9 chambers;
* content(11) = $2^{{11}}\\,5^3\\,11^2$ — support ⊆ primes$(n(n-1)) = \\{{2,5,11\\}}$ ✓, with $5^3$ ($5\\,|\\,10$)
  and $11^2$ ($11\\,|\\,11$) both honoring the refined patterns. Nine data points; the chain stands.

---

## 2. Cusps — the whisker tip

$\\gcd(p_{{10}}', p_{{10}}'') = 1$ exactly ⟹ **9 simple ordinary cusps** = roots of $p_{{10}}'$,
**EXACTLY 1 real** (Sturm — the reality dance $1,2,1,2,1,2,1,2,1$ extended to $d=10$, all exact).
The real cusp sits at $t = 0.3303841766\\ldots$ → $(s, r) = (0.321272,\\ 0.035372)$, and it is **NOT HIT**:
the residual octic has **0/8 real roots** (110-digit synthetic-division audits ≤ 4.2e-109) — so the
missed whisker ends in a cusp swallowtail tip, exactly the odd-fiber parity mechanism.

## 3. Nodes — the dance says four acnodes, zero crunodes

Bitangent Gröbner eliminant: degree $(n-1)(n-2) = 90$, and

$$\\mathrm{{elim}} = (55\\,p_{{10}}')^2\\cdot\\mathrm{{cofactor}}[72],\\qquad \\boxed{{K = 55^2 = 3025}}\\text{{ exact}},$$

cofactor **squarefree** and **coprime to $p_{{10}}'$** (exact gcds over ℚ). 72 roots at 120 digits
(max residue 2.5e-118), pairing with margins $10^{{116}}$: **36 nodes**, no triple points, no node-cusp
overlaps, budget $90 = 18 + 72$ ✓.

* **crunodes: 0** — the dance $0,1,0,1,0,1 \\to 0$ for $n = 11$ (whisker chambers have no cone corner);
* **acnodes: 4** — the staircase $1,1,2,2,3,3 \\to 4$:
  $(-3.5335, 2.4380)$, $(-0.2245, 1.7113)$, $(1.0265, 0.1411)$, $(-0.8422, -1.1019)$,
  each with 1/7 real residual roots (odd septic — the acnode modality).

---

## 4. The real side — whiskers don't miss open sets

* **census (200k, $\\sigma = 1.5$)**: {{1: **82.647%**, 3: **17.354%**, 5: **0%**}} — inside the locks
  [80.5, 84.5] / [15.5, 19.5] / [0, 0.5]. Compare the cone chamber d = 9: {{0: 8.726, 2: 84.03, 4: 7.25}}.
  The parity census law marches on.
* **200×200 region grid**: counts {{1, 3}} only — no empty open set exists; panel (a) shows the wall
  (crimson) crawling exactly along the count boundary to the whisker tip.
* **fiber counts**: generic 11 / fold 9 / cusp 8 / acnode **7** bounded — the node double-contact
  escapes twice (γ = 0 at both contacts), as budgeted.
* **escape laws**: fold $|\\gamma| \\sim 0.557\\,\\delta^{{0.4989}}$, cusp $\\sim 1.28\\,\\delta^{{0.6638}}$ —
  exponents $1/2$ and $2/3$, again universal ($A_2$, $A_3$).
* **hinge**: $C \\to 0$ antipodal pair **±0.90603** — the tower identity $p_1 = 2q_2$
  ($\\tfrac{{107}}{{55}} = 2\\cdot\\tfrac{{107}}{{110}}$) re-verified on its own terms.

---

## 5. Monodromy — $S_{{11}}$, Jordan style

Four loops, refinement-agreeing permutations of the 11 sheets, $\\min|D_{{10}}| \\ge 7.8e{{30}}$ on all loops:

| loop | permutation |
|---|---|
| fold | $(9,10)$ — **transposition** |
| cusp | $(9,11,10)$ — 3-cycle |
| $s = 200e^{{it}}$ | $(1,3,4,5,6,8,10,9,7,2)$ — 10-cycle |
| $r = 200e^{{it}}$ | $(1,2,4,5,7,8,10,11,9,6,3)$ — 11-cycle |

Transitive ✓, 2-homogeneous (**55/55** pairs reached by the verified generators) ⟹ primitive;
primitive + transposition ⟹ **$G = S_{{11}}$ by C. Jordan's theorem**. Bonus partial closure
$|G| \\ge 400{{,}}006$. The monodromy column of the tower saga now reads $S_3, S_4, S_5, S_6, S_7,
S_8, S_9, S_{{10}}, S_{{11}}$ — every chamber full-symmetric.

---

## 6. The figure

![The Whisker Chamber](data:image/png;base64,{fig})

**(a)** whisker atlas with wall, cusp tip and the four acnodes; **(b)** term law with the exact support
map (4 holes marked — the fourth fictitious — and $s^{{11}}$ out of the cone); **(c)** parity censuses
side by side and the eliminant's arithmetic; **(d)** monodromy + scoreboard.

---

## 7. HONESTY LEDGER

1. **Support-set cosmetic bug.** First stage-A holes print compared monomial tuples to stored
   *`(tuple, coeff)`* pairs and reported the whole 67-point cone as "missing" — the wall itself, terms
   count 64, and all certificates were unaffected. Re-derived directly: holes = expected three plus
   the (then-new) r¹¹-fictitious point, which *became* tonight's refinement of the term law. A bug
   that paid for itself.
2. **Self-ported Sturm routine undercounted** (dance printed 1,1,0,1,0,… against the known exact
   1,2,1,2,…) — replaced verbatim with note 10's proven implementation; all 9 counts then agreed
   with both Sturm and 110-digit `nroots`.
3. **Class-check transcription dropped the $-w\\,p(w)$ term** of the tangent-through-point equation —
   claimed lc should be $-(n-1)/n$; the true leading coefficient is $(1-n)\\,\\mathrm{{lead}}(\\Phi) = 10/11$.
   Fixed; certificate True.
4. **`complex()` vs mpmath string form** "(x + 0.0j)" — one more parse helper, two more crash cycles.
5. House rules honored elsewhere: no raw `polyroots` near multiple roots (synthetic deflation with
   remainder audits), gcds tested by degree, resultant-derived cofactors primitivized before
   squarefree verdicts.

---

## 8. Scoreboard

| # | Prediction (locked before computation) | Verdict |
|---|---|---|
| 1 | terms(11) = 64; holes {{const, s¹, r⁹}} | 🟢 exact support law on 2 walls |
| 2 | D10 degree 11, irreducible, param identity | 🟢 |
| 3 | D10(0,0) = 0; s² ∣ D10(s,0) | 🟢 |
| 4 | magnitude law L = 110, ratio 1 | 🟢 |
| 5 | content 2¹¹5³11²; support ⊆ {{2,5,11}}; 11² (11∣n), 5³ (5∣n−1) | 🟢 |
| 6 | fingerprint sign (−1)¹¹; magnitude = 2²⁰5²⁸11⁹ | 🟢 9/9 |
| 7 | 9 cusps; Sturm exactly 1 real (dance) | 🟢 |
| 8 | cusp NOT HIT — residual octic 0 real (whisker) | 🟢 |
| 9 | eliminant = (55p′)²·cofactor[72], K = 3025 | 🟢 exact |
| 10 | cofactor squarefree + coprime | 🟢 |
| 11 | 36 nodes = 0 crunodes + 4 acnodes + 32 complex | 🟢 |
| 12 | budget (n−1)(n−2) = 90 = 18 + 72 | 🟢 |
| 13 | fiber counts 11/9/8/7 | 🟢 |
| 14 | escape slopes 1/2, 2/3 | 🟢 0.4989, 0.6638 |
| 15 | census {{1, 3, 5}} in [80.5,84.5]/[15.5,19.5]/[0,0.5] | 🟢 82.65/17.35/0 |
| 16 | hinge antipodal (p₁ = 2q₂) | 🟢 ±0.90603 |
| 17 | det JF = 1 at 5/5 (d = 10) | 🟢 |
| 18 | monodromy $S_{{11}}$ (Jordan) | 🟢 + |G| ≥ 400006 |

**18 / 18 green.**

---

## 9. Queue (the train rolls on)

1. **Chamber n = 12 (d = 11)** — the even-fiber return of the cone. Locks: terms(12) = 78 − 2 = 76;
   K = den$(p_{{11}}')^2 = 22^2 = 484$; cusps 10, Sturm 2 (both hit, parity); nodes 45 = 1 crunode
   + 4 acnodes; census-cone per note-13's envelope law $m(11) = 8.7556\\%$ ± MC noise (lock
   [8.67%, 8.84%]); corner ∈ (note-12 lock) $(-0.950, -0.937)$, now sharpenable by the tie-model;
   monodromy $S_{{12}}$.
2. **The corner $\\sigma(d)$ law** vs the $-0.9733$ fit (note 12): the tie-point model from note 13
   meets the crunode migration series — $t_2(d) \\to -1$ like $-1 - c'/d$?
3. **Second-order closure** of the $(\\ln d)^2/d^2$ coefficient (note 13's amber, promoted to theory).
4. General-n proofs: $C_2/C_3$ closed forms; content-chain law (9 points, prime-support law solid;
   2-adic/3-adic exponents open).
5. Far lock frozen: $m(251) = 8.634481\\%$.
6. Swallowtail (5,1)-contact 4/5-slope hunt; generic-seed deg-5 control atlas; Moh 2-D chamber.

*Eleven chambers in, the tower has stopped surprising us and started keeping its word: every law,
every hole, every dance step — and the whisker, true to its name, mewsses nothing.* 🐈‍⬛🚂
"""

open("jacobian_undecic_atlas.md", "w").write(md)
print(f"saved jacobian_undecic_atlas.md ({len(md)/1024:.0f} KB)")
