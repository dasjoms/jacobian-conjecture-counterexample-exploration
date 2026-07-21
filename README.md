# Jacobian Conjecture Counterexample: AI Sandbox Exploration

> **An autonomous AI agent's structural investigation of the July 2026 Jacobian Conjecture counterexample, producing 23 lab notes, ~130 Python/Singular scripts, and a complete atlas of the seed-family tower from fiber degree 3 through 12 — including a proven eliminant factorization theorem (the Budget Law), a transversality dictionary reducing MAX-SING to two exact gcd computations, and a combinatorial proof of the surjectivity threshold at fiber degree 5 (the Threshold Deck).**

---

## What This Is

This repository contains the complete workspace output of an AI agent (Claude Fable 5) that was given free rein to explore the Jacobian Conjecture counterexample announced by [Levent Alpöge on July 19, 2026](https://x.com/__alpoge__/status/2079028340955197566). The agent ran in a sandboxed environment over approximately 24 hours, producing 23 sequential lab notes documenting its findings, along with all supporting code, data, and figures.

**This is AI-generated exploratory mathematics.** The agent was provided minimal human guidance and encouraged to follow its own curiosity. The work should be treated as a rich collection of machine-certified computational experiments and structural observations — not as peer-reviewed mathematical research. Some results have complete proofs; others are backed by extensive verification at finite parameters but remain conjectural in full generality.

---

## Background and Context

On July 19, 2026, mathematician Levent Alpöge (Anthropic) announced that Claude Fable 5 had found an explicit counterexample to the 87-year-old Jacobian Conjecture:

$$F(x,y,z) = \big((1+xy)^3 z + y^2(1+xy)(4+3xy),\; y + 3x(1+xy)^2 z + 3xy^2(4+3xy),\; 2x - 3x^2y - x^3z\big)$$

with $\det(JF) = -2$ identically, yet $F(0,0,-\tfrac14) = F(1,-\tfrac32,\tfrac{13}{2}) = F(-1,\tfrac32,\tfrac{13}{2}) = (-\tfrac14, 0, 0)$.

This counterexample has been widely and independently verified. See:

- **Original announcement:** [Alpöge's X post](https://x.com/__alpoge__/status/2079028340955197566)
- **Explainer with seed family:** [jacobianfun.org — The Jacobian counterexample, explained](https://jacobianfun.org/jacobian-explained) (Alexis Gallagher)
- **Galois/monodromy analysis:** [MathOverflow — Galois structure of the new counterexample](https://mathoverflow.net/questions/513387) (Qiaochu Yuan)
- **Geometric construction:** [Secret Blogging Seminar — The new counterexample to the Jacobian conjecture](https://sbseminar.wordpress.com/2026/07/20/the-new-counterexample-to-the-jacobian-conjecture/) (Will Sawin, Jordan Ellenberg)
- **News coverage:** [New Scientist — AI's solution to 87-year-old riddle](https://www.newscientist.com/article/2580374-ais-solution-to-87-year-old-riddle-takes-mathematicians-by-surprise/)
- **Wikipedia:** [Jacobian conjecture](https://en.wikipedia.org/wiki/Jacobian_conjecture) (updated July 2026)
- **Lean formalization:** [Formal Conjectures PR](https://github.com/google-deepmind/formal-conjectures) by Paul Rouzeau (Imperial College London)
- **MathWorld:** [Jacobian Conjecture](https://mathworld.wolfram.com/JacobianConjecture.html) (updated July 2026)
- **Verification preprint:** [arXiv — A Counterexample to the Jacobian Conjecture](https://arxiv.org/) (Zhang, 2026)
- **Direct consequences:** [zzhang-iu.github.io — Direct Consequences of the Three-Dimensional Counterexample](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html)

**What is widely reported and independently verified as of July 21, 2026** *(this list reflects the state of public discourse at the time of this repository's creation and may be incomplete — the mathematical community is moving quickly)*:

- The map has component degrees (7, 6, 4) and $\det J = -2$ (exact polynomial identity)
- Three distinct rational points collide at $(-1/4, 0, 0)$
- The map is not proper; non-properness does not invalidate the counterexample (properness is not a hypothesis of the conjecture)
- The conjecture is false for $n \geq 3$; $n = 2$ remains open
- The seed family $p_d(w) = 2w - 3w^2 + w(1-w)(w^{d-2} - 6/(d(d+1)))$ generates counterexamples for every fiber degree $n \geq 3$ (jacobianfun.org)
- Monodromy of the original map is $S_3$ (MathOverflow)
- A geometric construction via $\mathbb{P}^1 \times \text{Sym}^2(\mathbb{P}^1) \to \text{Sym}^3(\mathbb{P}^1)$ exists (Secret Blogging Seminar)

---

## Key Results

The following are the most significant findings from the agent's 23 lab notes, organized by confidence level. The "Builds on / Relates to" column describes what each result connects to — either prior public work, classical theory, or other results in this corpus. **We make no claims of novelty**; the mathematical literature is vast, and our search of prior work was necessarily incomplete. If any result has a known predecessor we missed, please open an issue.

### Results with Complete Proofs

These results have rigorous proofs (not just machine certification at finite parameters):

| # | Result | Note(s) | Builds on / Relates to |
|---|--------|---------|----------------------|
| 1 | **The Canal Theorem:** Maps of the form $F(x,y) = (A(xy)/x^m, B(xy)/x^n)$ have $\det JF = (nA'B - mAB')/x^{m+n}$, which can never be a nonzero constant for any $A, B$ (polynomial, entire, or formal). The 2D obstruction for the weighted-lift recipe is the single-channel structure, not polynomiality. | [Note 18](jacobian_missing_braid.md) | The 2D question was raised in public discourse but not answered; the recipe structure comes from the seed family (jacobianfun.org) |
| 2 | **The Shadow Uniqueness Theorem:** The diagonal $r = s$ meets the wall at exactly 3 points (odd $d$) or 2 points (even $d$), via the key identity $G'_d(t) = -(t-1)p'_d(t)$. The shadow root is unique; the pin's contact order is exactly 2. | [Note 17](jacobian_diagonal_census.md) | Builds on the shadow identity from the corner geometry analysis (Note 15) |
| 3 | **The Budget Law (divisibility):** The bitangent eliminant factors as $E_d = {p'_d}^2 \cdot \text{Cof}_d$ over $\mathbb{Q}$ for all $d \geq 2$, via the diagonal lemma ($(t_2-t_1)^2 \mid N$) and resultant multiplicativity. The cusp toll: every ordinary wall cusp pays exactly multiplicity 2 in the eliminant, with three transverse branches (slopes $-1, +1, -2$) at every flex. | [Note 22](jacobian_budget_law.md) | Builds on the eliminant column observed empirically in Notes 10–20; uses classical resultant theory and the Poisson formula |
| 4 | **Gal = $S_n$ for $d = 2..12$:** Three exact certificates (C1: transposition at the pin, C2: $(n-1)$-cycle at $s = \infty$, C3: $n$-cycle at $r = \infty$) plus a group lemma prove full symmetric monodromy. | [Note 19](jacobian_pin_transposition.md) | Extends $S_3$ monodromy established on [MathOverflow](https://mathoverflow.net/questions/513387) (Qiaochu Yuan) |
| 5 | **The Un-Rescue Theorem:** Every normalized degree-4 seed has at least one unrescued real cusp; hence no surjective-over-$\mathbb{R}$ Keller lift exists in this family. | [Note 5](jacobian_realghost.md) | Builds on the seed family's real behavior; relates to the Pinchuk counterexample for the real Jacobian conjecture |
| 6 | **The Connected-Cover Lemma:** A surjective, proper Keller map must be an automorphism; hence surjective non-injective Keller maps are forced to fail properness. | [Note 6](jacobian_escape_atlas.md) | Standard covering space theory applied to this setting; relates to Jelonek's non-properness framework |
| 7 | **The Pin Theorem:** $(-1, -1)$ lies on every wall; the diagonal is tangent at the pin (slope 1). | [Notes 15, 16](jacobian_corner_address.md) | Standalone structural observation about the wall family |
| 8 | **The No-Even-Shadow Theorem:** Even chambers provably have no near-corner diagonal crossing anywhere on $\mathbb{R}^-$. | [Note 16](jacobian_promotion.md) | Builds on the pin theorem and shadow identity |
| 9 | **The Surjectivity Threshold:** Maps with fiber degree $\geq 5$ in the family are surjective over $\mathbb{C}$, conditional on SQFREE ∧ T1 ∧ T2 (certified $d = 4..12$, chambers $n = 5..13$). The **partition argument** (Theorem A) proves this is combinatorially inevitable: every powerful fiber (all roots multiple) requires a partition of $n$ into parts $\geq 2$, and for $n \geq 5$ every such partition is killed by one of the three certificates — SQFREE kills parts $\geq 4$, T2 kills flex+contact overlaps, T1 kills tritangents. | [Notes 4, 23](jacobian_threshold_deck.md) | The partition argument is a standalone combinatorial theorem; the three certificates are the Budget Law's transversality (result #14 below); relates to Jelonek's non-properness theory |
| 10 | **The Two Expired Chambers:** The only non-surjective chambers are exactly identified: $n = 3$ ($d = 2$) has a unique perfect-cube fiber $h(w; 1/3, 1/27) = -(w-1/3)^3$ missing one rational curve; $n = 4$ ($d = 3$, Alpöge's chamber) has a unique perfect-square fiber $h(w; -1, -1) = -\frac{1}{4}(w^2+w-2)^2$ missing one rational curve. Gröbner basis confirms empty fiber over $(-1,-1,1)$. | [Note 23](jacobian_threshold_deck.md) | Exact unconditional results; the $n = 4$ curve is the missing curve from Notes 2–3; relates to the Budget Law's eliminant structure |
| 11 | **The Flat-Sheet Uniform Lemma:** The restriction to $\{x = 0\}$ is an elementary triangular automorphism for the entire recipe family, uniformly. The closed form $\lambda_d = m/(3(m-2))$ where $m = d(d+1)$ is proven for all $d$. | [Notes 19, 23](jacobian_pin_transposition.md) | The $C = 0$ fiber structure was partially analyzed in earlier notes; this unifies it |
| 12 | **The Hinge Antipodal Theorem:** The $C = 0$ gamma-sheet preimages are always antipodal in $x$, for every seed, every chamber. | [Note 5](jacobian_realghost.md) | Standalone identity about the $C = 0$ fiber structure |

### Results with Extensive Machine Certification (Conjectural in Full Generality)

These are backed by exact verification at many finite values of $d$ but not proven for all $d$:

| # | Result | Note(s) | Builds on / Relates to |
|---|--------|---------|----------------------|
| 13 | **The Max-Singularity Pattern:** Every wall has exactly $(n-2)$ cusps + $(n-2)(n-3)/2$ nodes = $(n-1)(n-2)/2$ delta-budget. Verified for $n = 3..12$. Equivalently: the contact system is transverse. | [Notes 6–11, 14, 20](jacobian_sextic_atlas.md) | The delta-budget framework is classical (genus formula, Plücker relations); the pattern that all singularities are ordinary is the empirical observation |
| 14 | **The Budget Law (transversality half):** T1 (no tritangents) ⟺ $\text{Cof}_d$ squarefree; T2 (no node-cusp overlap) ⟺ $\gcd(\text{Cof}_d, p'_d) = 1$; both ⟺ MAX-SING. Certified exactly for $d = 4..12$ (chambers $n = 5..13$) via two exact rational gcds. The all-$d$ transversality proof remains open, with three identified routes. | [Note 22](jacobian_budget_law.md) | The divisibility half is proven for all $d$ (result #3 above); this is the certified companion. The dictionary collapses numerical collision hunts into exact gcd computations |
| 15 | **The Envelope/Dual Curve Theorem:** The wall is the dual curve of $w \mapsto (w, \Phi(w))$; cusps = flexes, nodes = bitangents. Class = $n$ for every seed. | [Note 10](jacobian_nonic_atlas.md) | Classical dual curve / Legendre transform theory applied to this specific family |
| 16 | **The Escape Law $(m-1)/m$:** At a multiplicity-$m$ contact, preimages escape at rate $\delta^{-(m-1)/m}$. Measured at $m = 2, 3, 4$. | [Notes 6, 8, 10, 11](jacobian_generic_atlas.md) | Derives from singularity theory ($A_2$, $A_3$ normal forms); the universality across chambers is the empirical observation |
| 17 | **The Term Law:** $\text{terms}(n) = n(n+1)/2 - 1 - [n \geq 7]$, with the third hole $r^{n-2}$ opening at $n = 7$ via the $\beta$-theorem. | [Notes 10, 11](jacobian_nonic_atlas.md) | Builds on the discriminant support structure; the weighted-homogeneity of the discriminant is classical |
| 18 | **The Magnitude Law:** $\|[s^n] \text{resultant}\| = L^{2n-1}(n-1)^{n-1}\|a_n\|^{n-1}$, ratio exactly 1 in all chambers. | [Notes 10, 11](jacobian_nonic_atlas.md) | Relates to classical discriminant leading-term formulas |
| 19 | **The Eliminant Leading-Coefficient Laws:** $\text{lc}(E_d) = (-1)^d(d/(d+1))^{d-1}$, $\text{lc}(\text{Cof}_d) = (-1)^d d^{d-3}/(d+1)^{d-1}$, $\text{den}(p'_d) = d(d+1)/\gcd(d(d+1),6)$ — the last proven for all $d$. | [Note 22](jacobian_budget_law.md) | Builds on the magnitude law and content chain; the den law has a one-line proof, the lc laws are certified for $d = 4..13$ |
| 20 | **The Limit Shadow:** Tower seeds converge coefficient-wise to $2w - 3w^2$; walls converge to the original cubic's wall; missed-cone mass converges to $m^* = 8.6444647\%$ with boundary-layer corrections. | [Notes 12, 13](jacobian_limit_shadow.md) | Standalone asymptotic analysis of the seed family |
| 21 | **The Ringing Constant:** $m(d)$ rebounds toward $m^*$ with correction $\sim [A_0 - B(\ln d - \ln\ln d)]/d$. | [Note 13](jacobian_ringing.md) | Builds on the limit shadow; the logarithmic correction structure is the empirical observation |
| 22 | **The Reality Dance:** Real cusps per chamber: $1, 2, 1, 2, 1, 2, 1, 2, 1$ (Sturm-exact for $d = 2..10$). Crunodes: $0, 1, 0, 1, 0, 1, 0, 1, 0$. | [Notes 7–11, 14, 20](jacobian_septic_atlas.md) | Sturm theory applied to the seed family's $p'_d$; the alternating pattern is the empirical observation |
| 23 | **The Content Chain:** Prime support of $\text{content}(n) \subseteq \text{primes}(n(n-1))$; big primes $p \geq 5$: exponent 2 if $p \mid n$, 3 if $p \mid n-1$. | [Notes 10, 11](jacobian_decic_atlas.md) | Standalone arithmetic observation about the resultant's content |
| 24 | **The Parity Census:** Even $n$: missed open cones (~8.7%); odd $n$: whiskers only (measure zero). | [Notes 6–11, 14, 20](jacobian_octic_atlas.md) | Builds on the real behavior analysis; relates to the fiber-degree parity mechanism |

### Results That Extend or Connect to Known Public Work

| # | Result | Note(s) | Public Source |
|---|--------|---------|---------------|
| 25 | **Seed family counterexamples for all $n \geq 3$** (with geometry, not just existence) | [Note 1](jacobian_lab_notes.md) | [jacobianfun.org](https://jacobianfun.org/jacobian-explained) identified the family and listed counterexamples for degrees 3–100 |
| 26 | **$S_3$ monodromy of the original map** | [Note 4](jacobian_surjectivity.md) | [MathOverflow](https://mathoverflow.net/questions/513387) (Qiaochu Yuan) |
| 28 | **Basic verification (det, collision, escape curve)** | [Notes 1–3](jacobian_lab_notes.md) | Widely and independently verified by the community |

---

## Lab Notes Guide

The 23 lab notes are sequential. Each builds on the previous. The agent maintained "honesty ledgers" in every note documenting caught bugs, corrected errors, and falsified predictions — **please read these**. They are the agent's own quality control mechanism.

| # | File | Title | Key Content |
|---|------|-------|-------------|
| 1 | `jacobian_lab_notes.md` | The Jacobian counterexample — sandbox lab notes | Initial verification, factored form, fiber equation, homemade counterexample H |
| 2 | `jacobian_anatomy.md` | Anatomy of a counterexample | Missing curve $M$, wall hypersurface, stratification (3/1/0 fibers) |
| 3 | `jacobian_verification.md` | Trust, but verify | Three-method cross-check: structure, $\mathbb{F}_{101}^3$ census, Singular Gröbner cover |
| 4 | `jacobian_surjectivity.md` | Monodromy, missing curves, and surjective non-invertible Keller maps | $S_3$ monodromy, surjectivity threshold at fiber degree 5, dimension counting |
| 5 | `jacobian_realghost.md` | The real ghost: cusp obstructions and the un-rescue theorem | Real behavior, un-rescue theorem, Nagata conjugation, finite-field census |
| 6 | `jacobian_escape_atlas.md` | The escape atlas: wall, catastrophes, and full braid of $F_4$ | Wall of $F_4$ (quintic, 3+3=6), $S_5$ monodromy, escape rates |
| 7 | `jacobian_sextic_atlas.md` | The sextic chamber | Fiber 6: 4 cusps + 6 nodes = 10, $S_6$, cone 8.69%, antipodal hinge |
| 8 | `jacobian_septic_atlas.md` | The septic chamber | Fiber 7: 5 cusps + 10 nodes = 15, $S_7$, whisker, escape universality |
| 9 | `jacobian_octic_atlas.md` | The octic chamber | Fiber 8: 6 cusps + 15 nodes = 21, $S_8$, cone 8.74%, node-whisker |
| 10 | `jacobian_nonic_atlas.md` | The nonic chamber — and the Plücker key | Fiber 9: 7+21=28, $S_9$, envelope theorem, term law, magnitude law |
| 11 | `jacobian_decic_atlas.md` | The decic chamber | Fiber 10: 8+28=36, $S_{10}$, three-hole law, content chain, $A_3$ universality |
| 12 | `jacobian_limit_shadow.md` | The limit shadow | Tower convergence to original cubic, $m^* = 8.6444647\%$, boundary layers |
| 13 | `jacobian_ringing.md` | The ringing constant | $m(d)$ rebound, correction law, 20 exact data points to $d = 201$ |
| 14 | `jacobian_undecic_atlas.md` | The whisker chamber (fiber 11) | 9+36=45, $S_{11}$, 64 terms, all tower laws simultaneously verified |
| 15 | `jacobian_corner_address.md` | The corner's address | Pin theorem, shadow identity, asymptotic series, trio coalescence |
| 16 | `jacobian_promotion.md` | The promotion | Five numerical laws promoted to theorems: pin, dance, ghost, left cusp, shadow |
| 17 | `jacobian_diagonal_census.md` | The diagonal census | Shadow uniqueness theorem for all $d$, pin contact order = 2 |
| 18 | `jacobian_missing_braid.md` | The missing braid | **Canal Theorem**, 2D obstruction, braid ledger, entire avatar |
| 19 | `jacobian_pin_transposition.md` | The pin's transposition | Gal = $S_n$ proof ($d \leq 12$), flat-sheet lemma, errata from external reviewer |
| 20 | `jacobian_last_chamber.md` | The last chamber (fiber 12) | Completed atlas $n = 5..12$, eliminant column, corner-chases-shadow |
| 21 | `jacobian_generic_atlas.md` | The generic atlas | 8/8 random seeds, swallowtail singularity, escape law at $m = 4$ |
| 22 | `jacobian_budget_law.md` | **The budget law** | **Eliminant factorization theorem** (divisibility for all $d$, transversality dictionary), three lc/den laws, exact column through $d = 13$ |
| 23 | `jacobian_threshold_deck.md` | **The threshold deck** | **Surjectivity flagship consolidation**: partition argument (Theorem A), two expired chambers (Theorem B), flat sheet $\lambda$-law, literature positioning |

---

## Repository Structure

The workspace is **flat** — all 284 files live in one directory. Scripts reference data files by bare filename and must be run from the workspace root. The organizational structure is implicit in the filename prefixes:

**Lab Notes (23 markdown files)**
- `jacobian_lab_notes.md` through `jacobian_threshold_deck.md`
- Sequential narrative of the agent's investigation
- Each note has an honesty ledger documenting caught bugs and corrected errors
- See the [Lab Notes Guide](#lab-notes-guide) above for a summary of each

**Core Verification Scripts**
- `verify_jacobian.py` — Basic det + collision check (start here)
- `jacobian_structure.py` — Factored form, symmetry, escape curve
- `jacobian_anatomy*.py` — Missing curve, wall analysis (3 scripts)
- `jacobian_family_and_fibers.py` — Fiber equation framework
- `jacobian_homemade.py` — Homemade counterexample H
- `family_maps_verify.py` — Seed family verification

**Atlas Scripts (per chamber)**
- `jacobian_atlas{1,2,3,3b}.py` — Chambers n=3,4
- `jacobian_atlas5_*.py` — Chamber n=5 (7 scripts)
- `jacobian_atlas6_*.py` — Chamber n=6 (5 scripts)
- `jacobian_atlas7_*.py` — Chamber n=7 (5 scripts)
- `jacobian_atlas8_*.py` — Chamber n=8 (6 scripts)
- `jacobian_atlas9_*.py` — Chamber n=9 (8 scripts)
- `jacobian_atlas10_*.py` — Chamber n=10 (6 scripts)
- `jacobian_atlas11_*.py` — Chamber n=11 (5 scripts)
- Each chamber has: wall computation, node/cusp analysis, monodromy, properness, figure

**Budget Law Scripts**
- `jacobian_budgetlaw_1.py` — Universal algebra certificates (stage 1)
- `jacobian_budgetlaw_2.py` — Exact column + pairing exhibit (stage 2)
- `jacobian_budgetlaw_fig.py` — Figure generation
- Data: `budgetlaw_stage{1,2}.json`, figure `budgetlaw_figure.png`

**Flagship / Threshold Deck Scripts**
- `jacobian_flagship_1.py` — Surjectivity threshold verification
- `jacobian_flagship_fig.py` — Figure generation
- Data: `flagship_stage1.json`, figure `flagship_figure.png`

**Structural Analysis Scripts**
- `jacobian_monodromy.py` — S₃ monodromy computation
- `jacobian_census*.py` — Fiber census scripts
- `jacobian_modp_*.py` — Finite-field computations
- `jacobian_nagata.py` — Nagata automorphism conjugation
- `jacobian_realghost*.py` — Real behavior analysis (7 scripts)
- `jacobian_generic*.py` — Generic seed atlas (3 scripts)
- `jacobian_flatsheet.py` — Flat-sheet uniform lemma

**Asymptotic / Limit Scripts**
- `jacobian_limits_*.py` — Limit shadow (4 scripts)
- `jacobian_limitshadow_build.py` — Envelope machinery
- `jacobian_ring_*.py` — Ringing constant (7 scripts)
- `jacobian_corner_*.py` — Corner geometry (8 scripts)

**Theorem Proving Scripts**
- `jacobian_promote_*.py` — Promotion to theorems (4 scripts)
- `jacobian_diagcensus_*.py` — Diagonal census (2 scripts)
- `jacobian_galois_*.py` — Galois = Sₙ certificates (2 scripts)
- `jacobian_holes_0.py` — Term law / β-theorem

**2D Investigation Scripts**
- `jacobian_keller2d_*.py` — Canal Theorem (3 scripts)
- `jacobian_keller2d_fig.py` — 2D figure

**Singular Scripts (6 files)**
- `grobcov_F.sing` — Gröbner cover of F
- `grobcov_F2.sing` — Gröbner cover (variant)
- `missing_curve.sing` — Missing curve M computation
- `surjectivity.sing` — Surjectivity proof
- `surjectivity_setup.sing` — Surjectivity setup
- `verify_M.sing` — Verify missing curve

**Data Files (~71 JSON, ~27 text, ~11 NumPy)**
- `atlas{5..11}_*.json` — Per-chamber data (walls, nodes, cusps, escape rates, monodromy, census)
- `atlas{4..11}_wall.txt` — Wall equations (verbatim)
- `atlas{5..11}_elim_raw.txt` — Raw eliminant polynomials
- `atlas{8..11}_cofactor.txt` — Cofactor polynomials
- `atlas{8..11}_monodromy.txt` — Monodromy permutations
- `atlas{5..11}_bitangent_eliminant.txt` — Bitangent eliminants
- `atlas_eliminant_column.json` — Eliminant identity across all chambers
- `jcorner_*.json` — Corner geometry data
- `limits_*.json` — Limit shadow data
- `promote_stage*.json` — Promotion certificates
- `ring_*.json` — Ringing constant data
- `diagcensus_stage1.json` — Diagonal census certificates
- `flatsheet_stage.json` — Flat-sheet lemma data
- `galois_exact_stage.json` — Galois certificates
- `keller2d_stage*.json` — 2D Canal Theorem data
- `generic_atlas.json` — Generic seed data
- `holes0_beta.json` — β-theorem data
- `swallowtail_*.json` — Swallowtail singularity data
- `atlas{5..11}_grid.npz` — Real atlas grid data
- `jcorner_map_d{2..5}.npz` — Corner map evaluations

**Figures (22 PNG files)**
- `anatomy_{F,G}.png` — Fiber anatomy
- `atlas{8,10,11}_figure.png` — Chamber atlases
- `atlas8_atlas.png` — Octic atlas
- `budgetlaw_figure.png` — Budget law (incidence geometry)
- `flagship_figure.png` — Threshold deck (surjectivity)
- `diagcensus_figure.png` — Diagonal census
- `escape_atlas.png` — Escape geometry
- `galois_figure.png` — Galois certificates
- `generic_atlas.png` — Generic seeds
- `jacobian_limits_figure.png` — Limit shadow
- `jcorner_figure.png` — Corner geometry
- `keller2d_figure.png` — 2D Canal Theorem
- `monodromy.png` — Monodromy visualization
- `{octic,septic,sextic}_atlas.png` — Chamber walls
- `promote_figure.png` — Promotion theorems
- `realghost.png` — Real behavior
- `ring_figure.png` — Ringing constant

---
## How This Agent Was Steered

**The honest answer: barely at all.** This workspace is overwhelmingly the product of the model's own "scientific curiosity." Here is the complete record of human intervention:

1. **Initial prompt (Notes 1–3):** The agent was given the counterexample and asked to verify it. From there, it autonomously began exploring the map's structure, discovering the factored normal form, odd symmetry, escape curve, fiber equation, and missing curve.

2. **Free exploration (Notes 4–15):** After the initial verification, the agent was told to freely experiment with whatever it found most interesting. It chose to:
   - Build a complete atlas of the seed family tower (fiber degrees 5 through 11)
   - Discover the max-singularity pattern, escape law, and monodromy tower
   - Investigate real behavior, the limit shadow, and asymptotic constants
   - Prove the un-rescue theorem and surjectivity threshold
   - Study the corner geometry and shadow identity

   During this phase, the agent began "climbing the tower" — computing one new chamber per turn.

3. **The 2D nudge (Note 18):** After Note 17, the human operator nudged the agent toward the missing 2D question — the one chamber it hadn't explored. The agent logged this nudge explicitly in Note 18 (`jacobian_missing_braid.md`): *"The nudge, and the chamber... the gracious conductor of this lab."* The result was the Canal Theorem.

4. **External review (Notes 19–23):** After Note 17, a snapshot of the workspace was handed to a **separate instance of Claude Fable 5** in an independent sandbox and tasked with:
   - Evaluating the original agent's findings
   - Verifying claims using its own code (not the original agent's scripts)
   - Proposing topics for further investigation

   The reviewer's notes were handed back to the original agent, which incorporated them starting with Note 19 (`jacobian_pin_transposition.md`). The agent explicitly credits the "external reviewer" throughout Notes 19–23.

5. **No other intervention.** The human operator did not suggest specific theorems, correct errors, or direct the mathematical content. The agent's honesty ledgers document its own caught bugs, false starts, and corrected mistakes.

---

## Technical Notes

**Software environment:** The agent used Python 3 with SymPy 1.14 (exact symbolic computation), mpmath (50–200 digit arithmetic), NumPy, and Matplotlib. Singular 4.4.0 was used for Gröbner basis computations (Notes 3, 4).

**Verification methodology:** The agent employed multiple independent verification strategies:
- Exact symbolic computation over ℚ (SymPy)
- Exhaustive finite-field enumeration over $\mathbb{F}_{101}^3$ (1,030,301 points)
- Independent CAS verification (Singular)
- High-precision numerical audits (50–200 digits via mpmath)
- Locked predictions registered in script headers before computation

**Known limitations:**
- Some "theorems" are backed by machine certification at finite $d$ (e.g., $d \leq 300$) plus induction arguments, not traditional all-$d$ proofs
- The asymptotic series for the shadow identity may be divergent ($|u_5| \approx 1$)
- The agent caught and documented its own errors in honesty ledgers — read these
- The external reviewer's input has only been incorporated in Notes 19 and 20

---

## Inviting Collaboration

This work is released in the hope that the mathematical community will:

1. **Verify** the claims independently — especially the Canal Theorem, the max-singularity pattern, and the Galois = $S_n$ proof
2. **Extend** results to all $d$ — the corridors to all-$d$ proofs are identified but not completed
3. **Formalize** the strongest theorems in Lean, Coq, or similar systems
4. **Connect** the seed family to geometric constructions (Sawin's $\mathbb{P}^1 \times \text{Sym}^2(\mathbb{P}^1)$ picture)
5. **Explore** the open questions listed in the lab notes: the Budget Law's transversality corridor (three identified routes: swap-symmetry induction, mod-$p$ certified spread, the B10 coefficient ladder), the all-$d$ Galois proof, the content chain, and Moh's 2D chamber

If you find errors, please open an issue. If you extend the work, please share. If you know of prior work that any of these results overlap with, please let us know so we can credit it properly.

---

## License

This work is released under the [MIT License](LICENSE). The mathematical content is factual and may be freely used, verified, and extended.

---

*Generated by a Claude Fable 5 agent running autonomously in a sandboxed environment, July 20–21, 2026. The agent was provided the Alpöge counterexample and minimal guidance, then encouraged to follow its curiosity. 284 files, 23 lab notes, ~130 scripts, 22 figures.*
