# certified-packing

Protein side-chain packing with a **machine-checkable proof of optimality**.

Side-chain packing — choose one rotamer per residue on a fixed backbone to
minimise a pairwise energy — is a discrete pairwise MRF. This repository
certifies its global optimum with a certificate a third party can verify **in
exact rational arithmetic, without rerunning any solver**.

There are two certificate routes, and the second supersedes the first:
- **The QUBO route.** Pose the MRF as a QUBO and solve it with
  [CROWN](https://github.com/codenlighten/crown). It certified 30 of 42
  benchmark instances, to a floating-point tolerance, and fails when the
  residual core's treewidth exceeds ~27.
- **The LP-dual route.** The local-polytope LP plus triplet clusters, solved
  exactly by SoPlex and checked with Python `Fraction`s. It certifies **all
  42** exactly (gap 0) on the DEE-reduced instances, and **38 of them** on the
  unpruned instances, with no DEE in the chain. See
  [Exact certificates for all 42](#exact-certificates-for-all-42-the-relaxation-not-the-core).

The resulting statement is categorically different from what annealing gives:

> **"This packing is the proven global minimum of the stated energy function"**
> — not "this is the best my sampler found."

**Read this before the tables below.** The claim is about the *artifact*, not
about speed or scale. Benchmarked head-to-head, `toulbar2` solves every instance
here in under a second, including the twelve the CROWN route could not certify.
The exact LP certificates take seconds to minutes each. If you want an optimal
packing, use `toulbar2`. What this repository offers
is a certificate you can re-check without trusting the solver — see
[Measured against `toulbar2`](#measured-against-toulbar2).

## The result

The naive reduction fails at useful sizes. Goldstein **Dead-End Elimination**
fixes it, and the reason is structural rather than incidental. **The finding
holds under both a coarse counting score and a physics-based
molecular-mechanics energy** — see [Energy function](#energy-function).

| arm | certified | mean compression | mean vars |
|---|---|---|---|
| A — one-hot (baseline) | **3/5** | 0.000 | 74 |
| B — DEE + one-hot | **5/5** | 0.656 | 25 |
| C — DEE + arity-aware | **5/5** | 0.377 | 7 |

Per instance, 4 rotamers per flexible residue:

| protein | residues | rotamers | DEE pruned (steric / MM) | A | B | C |
|---|---|---|---|---|---|---|
| 2I9M | 9 | 36 | 69% / 72% | ✔ | ✔ | ✔ |
| 1E0Q | 16 | 64 | 67% / 47% | ✔ | ✔ | ✔ |
| 1L2Y | 13 | 52 | 60% / 65% | ✔ | ✔ | ✔ |
| 1FME | 25 | 100 | 68% / **75%** | ✘ | ✔ | ✔ |
| 1VII | 30 | 120 | 68% / 72% | ✘ | ✔ | ✔ |

Arm results are identical under both energies: A certifies 3/5, B and C certify
5/5.

**On 1FME the uncertified baseline was 157.90 worse than the true optimum**
(678.40 vs 520.50) — and without a certificate there was no way to know. That is
the practical case for certification in one line.

## Why the encoding decides everything

- **Direct binary** (2 rotamers, one variable per residue, no penalty term):
  roof duality collapses the entire problem — compression 1.000, empty core.
- **One-hot** (k > 2, with `λ(Σx−1)²` per residue): roof duality fixes
  **nothing** — compression 0.000 in every case. The penalty places a `+2λ`
  coupling between every pair of rotamers within a residue: a dense, strongly
  non-submodular clique, precisely what roof duality cannot handle. Certification
  falls back to exact core solving, which is exponential in treewidth.

  On the original five small proteins this looked like a **variable-count**
  ceiling — it held to ~90 variables and failed at 100–120. The size ladder
  showed that reading was wrong: see
  [Larger proteins](#larger-proteins-the-ceiling-is-not-a-size-ceiling).
  Treewidth is the binding constraint, and variable count was only standing in
  for it on a sample too small to tell them apart.

**The fix is not a cleverer binary encoding. It is to apply persistency before a
binary encoding exists.** Goldstein DEE is categorical persistency — it prunes
rotamers that provably cannot appear in any global optimum, operating directly on
the MRF. It removes 60–69% of rotamers here, leaving most residues with one or
two survivors, at which point the penalty clique largely disappears.

DEE (Desmet et al. 1992) and roof duality are the same idea at different levels.
Running the categorical one first is what makes the binary one work.

## Energy function

Two energies are implemented (`packing/forcefield.py`), and the pipeline is run
against both so the conclusion cannot rest on an artefact of a toy score:

- **`steric`** — the coarse counting score: clashes minus contacts.
- **`mm`** — softened Lennard-Jones 6-12 with per-element radii and well depths
  (AMBER-like), a Gaussian hydrogen-bond term between polar atoms, and a
  burial-based solvation term. The LJ repulsion is linearised below 0.6·r_min,
  Rosetta-style, so clashes stay finite.

**The result is stable across both.** Under `mm`, DEE prunes 47–75% of rotamers
and certification still goes 3/5 → 5/5. Two details worth recording:

- **On 1FME, DEE alone solved the instance outright** — every residue reduced to
  a single surviving rotamer, leaving *zero* variables for the QUBO solver.
  Categorical persistency was sufficient on its own.
- **1E0Q is the hard case.** DEE pruned only 47%, leaving arities up to 4, and
  certification came from exact-core solving rather than roof duality
  (compression 0.147). It is the instance to watch as the energy gets more
  realistic.

Uncertified baselines under `mm` came within 10.76 and 4.33 kcal/mol of the
optimum — small in relative terms against total energies of ~4400, but 10 kcal/mol
is chemically significant, and without a certificate there is no way to know
which case you are in.

`mm` is a **simplified** molecular-mechanics energy: no electrostatics, no
torsional strain, no rotamer priors, and burial-count solvation rather than
EEF1/GB. It exists to test that the pipeline survives a continuous, physically
motivated energy with realistic dynamic range — not to predict structures.

## Rotamer arity — does it scale?

Published libraries carry 10–100 rotamers per residue. `packing/rotamers.py`
generates rotamers over multi-χ staggered wells (**not** the Dunbrack library,
which requires registration — see that module's docstring), driving mean arity
from ~4 to ~24. Under the `mm` energy:

| protein | χ≤1 mean k / pruned | χ≤2 mean k / pruned | χ≤3 mean k / pruned | certified |
|---|---|---|---|---|
| 2I9M (9 res) | 3.7 / 67% | 10.3 / 85% | 23.3 / **91%** | ✔ ✔ ✔ |
| 1L2Y (13) | 3.3 / 60% | 8.6 / 74% | 14.2 / **84%** | ✔ ✔ ✔ |
| 1E0Q (16) | 4.0 / 69% | 8.7 / 86% | 14.3 / **90%** | ✔ ✔ ✔ |
| 1FME (25) | 4.0 / 70% | 11.1 / 83% | 24.4 / **89%** | ✔ ✔ ✔ |
| 1VII (30) | 3.8 / 70% | 10.7 / 88% | 18.5 / **90%** | ✔ ✔ ✔ |

(Split DEE enabled; pruning percentages are post-split.)

**DEE gets stronger as arity grows, not weaker** — pruning rises from 60–70% at
χ≤1 to 81–90% at χ≤3. This is the opposite of the expected failure mode, and it
makes sense: more rotamers give each residue more chances to find a dominating
alternative, so Goldstein's criterion fires more often.

With Goldstein alone, 13/15 certified: 1FME and 1VII at χ≤3 pruned 85% and 81%
but still left **92 and 105 residual variables**, above the ~90-variable
certification ceiling.

*A prediction of ours that was wrong:* 1E0Q was flagged as the canary after
pruning only 47% at 4 rotamers. It scaled fine (69% → 81% → 85%). Pruning rate
on a small instance did not predict the limit; total residual size did.

### Split DEE closes the gap: 15/15

Split DEE (Pierce, Spriet, Desmet & Mayo 2000) partitions the environment by a
witness residue and allows a **different** dominating rotamer per partition —
strictly stronger than Goldstein, which needs one rotamer to dominate everywhere
at once. On the two failures:

| instance | Goldstein survivors | Split survivors | vars | result | time |
|---|---|---|---|---|---|
| 1FME χ≤3 | 92 | **66** (89% pruned) | 80 → **53** | bracket → **certified** | 61s → **10s** |
| 1VII χ≤3 | 105 | **56** (90% pruned) | 87 → **31** | bracket → **certified** | 130s → **7s** |

**All 15 instances now certify**, and the hard cases got *faster*: a smaller
residual core means the solver finishes quickly rather than grinding to a
bracket. Split DEE is verified optimum-preserving against brute force in
`tests/test_packing.py`.

## Larger proteins: the ceiling is not a size ceiling

Every instance in the sections above certified in under 10 seconds, which means
the size limit was never actually probed — those proteins were 9–30 flexible
residues. `experiments/size_scaling.py` runs the identical pipeline over a
ladder of 14 proteins, **32–369 flexible residues**, at all three χ depths.

**Certification reaches 369 flexible residues** — 12× the largest instance in
the original set — and 41 of 42 instances were small enough to attempt.

| χ depth | certified | largest certified | failures |
|---|---|---|---|
| χ≤1 | **14/14** | 1GAI, 369 residues | — |
| χ≤2 | **12/14** | 1ADE, 336 residues / 251 vars | 1A6M, 1GAI |
| χ≤3 | **4/13** | 1AKI, 103 residues | 9 (1 more skipped over the 400-var cap) |

At χ≤1 every instance certifies by roof duality alone (`bound-tight`), in ≤7.5s
including energy construction, all the way to 369 residues. DEE pruning does not
decay with size — it holds at 67–73% at χ≤1 and 81–89% at χ≤3 regardless of
whether the protein has 32 residues or 369.

### The prediction that failed, and what replaced it

The earlier sections recorded a **~90-variable certification ceiling**. The
ladder falsifies it in both directions:

| instance | residual vars | result |
|---|---|---|
| 1ADE χ≤2 | **251** | ✔ certified, exact-core, 57s |
| 3PGK χ≤2 | 163 | ✔ certified, exact-core, 20s |
| 8ABP χ≤2 | 150 | ✔ certified, exact-core, 32s |
| 1QOP χ≤2 | 94 | ✔ certified, **bound-tight, 1.4s** |
| 1UBQ χ≤3 | **93** | ✘ bracket |
| 1A6M χ≤2 | 133 | ✘ bracket |

A core of 251 variables certifies while one of 93 does not. Nor is it protein
size: **1AKI certifies at 103 flexible residues while 1UBQ fails at 65.**

Three candidate explanations were measured and all three fail, because the
certified and failed ranges overlap:

| quantity | certified | failed | verdict |
|---|---|---|---|
| flexible residues | 32–369 | 65–369 | no signal |
| residual variables | 2–251 | 93–354 | overlaps 93–251 |
| one-hot (arity ≥ 3) variables | 0–186 | 77–298 | overlaps 77–186 |

The "~90 variables" number was never a property of the method. It was an
artifact of five small proteins on which variable count happened to correlate
with the thing that actually binds.

### What actually binds: core treewidth

CROWN certifies in two stages. Roof duality fixes what it can; whatever survives
is a **core** that must be solved exactly, and exact solution is exponential in
the core's induced **treewidth** — not in its variable count. CROWN reports both
(`core_size`, `core_width`), and `size_scaling.py` now records them.

Treewidth separates where the other three quantities did not (41 attempted
instances, 30 certified):

| core treewidth | certified | failed |
|---|---|---|
| **< 22** | **27** | 0 |
| 22–27 | 3 | 2 |
| **> 27** | 0 | **9** |

Outside a six-unit transition band the split is total. Inside it, search
difficulty rather than width decides — CROWN's core solver is branch-and-bound
under a budget, not pure elimination, so width bounds the cost without fixing
it.

The decisive comparison, both instances real:

| instance | vars | core treewidth | result |
|---|---|---|---|
| 1ADE χ≤2 | **251** | 25 | ✔ certified, 53s |
| 1UBQ χ≤3 | **93** | 31 | ✘ bracket |

Nearly three times the variables, certified — because the core is narrower.
This is why **pruning harder is the wrong lever**: pair-split DEE (`split=2`)
removes more rotamers but leaves the interaction graph the same shape, and on
the two instances it was built for it moved 1UBQ from 93 to 87 variables and
3CHY from 98 to 96 without either crossing the line. What has to shrink is
width.

## Measured against `toulbar2`

The README has always conceded that exact solvers do this already. That
concession was never measured, which made it the load-bearing unknown in the
repository. `experiments/toulbar2_benchmark.py` measures it: the same MRF goes
to `toulbar2` — a weighted-CSP solver that takes the pairwise problem directly,
with no binary encoding — and to this pipeline. 42 instances, 32–369 flexible
residues, all three χ depths.

**The result is not close.**

| | toulbar2 | this pipeline |
|---|---|---|
| instances solved to optimality | **42 / 42** | 30 / 42 |
| slowest instance | **0.91 s** | 361 s |
| median time on instances we certify | 0.04 s | 0.9 s (**11× slower**) |
| worst case on instances we certify | — | **154× slower** |

`toulbar2` solved every instance in under a second, including all twelve this
pipeline cannot certify. On those twelve, our best answer was worse than the
true optimum by between 2.4 and **867.8** kcal/mol:

| instance | our shortfall | toulbar2 | ours |
|---|---|---|---|
| 1GAI χ≤3 | **+867.79** | 0.62 s | 361 s |
| 1A6M χ≤3 | +669.68 | 0.25 s | 174 s |
| 1QOP χ≤3 | +364.01 | 0.39 s | 107 s |
| 1UBQ χ≤3 | +57.21 | 0.09 s | 34 s |
| 1A6M χ≤2 | +2.42 | 0.09 s | 34 s |

### What the benchmark does confirm

**Every one of the 30 certified optima matches `toulbar2` exactly.** That is the
first independent check that these certificates are correct rather than merely
self-consistent — an encoding bug or an unsound DEE step would show up here as a
disagreement, and none did. It also confirms the pipeline never claimed
optimality it did not have: on all twelve failures it returned `bracket`, not a
false certificate. An annealer would have reported the 1GAI χ≤3 answer with no
indication it was 868 off.

### What it does not

It does not support any claim of speed or scale. On time-to-optimal `toulbar2`
dominates: faster on every non-trivial instance, and it solves a strictly larger
set. **If you want an optimal packing, use `toulbar2`, not this.**

What remains is narrower and worth stating precisely: this pipeline emits a
**certificate a third party can re-check by bounded arithmetic without rerunning
the solver**. `toulbar2` proves optimality internally — you trust the solver, or
you rerun it. That difference is the contribution, and it is an artifact
difference, not a performance one. A benchmark table cannot express it, which is
why the table above should be read as bounding the claim rather than supporting
it.

## Exact certificates for all 42: the relaxation, not the core

The treewidth section above says exact core solving is the wall. The way round it
is not a narrower core. It is to stop solving the core at all, and certify with
a **relaxation bound** instead: if a lower bound LB on every assignment's energy
meets the energy of a known assignment, that assignment is optimal, whatever the
treewidth.

The relaxation is the **local marginal polytope LP** of the categorical MRF:
one rotamer per residue as a hard constraint, not the one-hot penalty that
defeats roof duality. Its dual is a **reparameterisation**, a set of messages
that move energy between each pair table and its endpoints without changing any
assignment's total. For *any* messages, the sum of the reparameterised tables'
minima is a lower bound. So the certificate is an assignment plus the messages.
The **verifier** needs no solver: whatever produced the messages (SoPlex,
HiGHS, a bug) is untrusted input to it. A bad message can only weaken the bound,
never make it wrong.

| rung | what | result on the 42 |
|---|---|---|
| (a) | local-polytope LP | tight on **36/42** (to floating-point noise). The other 6 have real gaps of 0.001–0.97 kcal/mol |
| (b) | + triplet clusters (Sontag et al. 2008) on triangles touching fractional residues | closes **all 6** in one round, 24–196 triplets |
| (c) | + semidefinite constraints | **not needed** |

The six instances that needed triplets are 8ABP χ≤2, 1ADE χ≤2, 1GAI χ≤3,
1QOP χ≤3, 2LZM χ≤3 and 3PGK χ≤3.

### From "tight to 1e-12" to exact

A floating-point LP that is "tight to 1e-12" is not a proof. Three steps make it
one:

1. **Exact dual.** Each LP is written with exact rational coefficients (every
   stored double is a dyadic rational, written as `p/q`) and solved by
   [SoPlex](https://soplex.zib.de) 8.1.0 in exact rational mode (iterative
   refinement, Gleixner, Steffy & Wolter), which returns a rational dual.
   Nothing downstream trusts that dual; it only has to pass step 2.
2. **Exact check.** `packing/certify.py` re-derives LB from the messages with
   Python `Fraction`s only. It needs no LP solver, no floating point and no
   tolerance. It enumerates every entry of every table and calls the assignment
   optimal only if **LB ≥ E(x\*) exactly**. The assignment is SoPlex's exact
   primal solution, which was integral on all 42 instances and agreed with
   `toulbar2`'s on all 42. Its energy is recomputed exactly by the checker in
   any case.
3. **No DEE in the chain.** A certificate on the Split-DEE-reduced instance
   still trusts DEE, a floating-point procedure. `experiments/lift_certificate.py`
   lifts each reduced certificate to the **unpruned** instance:
   - Messages on kept rotamers are copied. For a pair (i, j) with reparameterised
     table P and minimum m, pruned rotamers get, in exact arithmetic:
     `d_j(b) = min over kept a of [P(a,b) − d_i(a)] − m`, then
     `d_i(a) = min over all b of [P(a,b) − d_j(b)] − m`.
     No entry of the pair term then falls below m. Triplet-message entries that
     touch a pruned rotamer are set low enough that no triplet minimum moves.
   - The few pruned rotamers that still fall short go back into a working set,
     the (still small) LP is re-solved exactly, and the lift repeats. This is
     row-and-column generation, with DEE's kept set only as the starting point.
   - The result is checked on the unpruned instance. **That passing check is
     the proof.** Neither the lifting rule nor DEE has to be correct for the
     claim to hold; a wrong lift would simply fail the check.

   Solving the unpruned LP exactly and directly works too: 1ADE χ≤1 and
   3PGK χ≤1 were done both ways, and both ways prove the same optimum. But it
   took 234–321 s at 40k columns, and the χ≤3 LPs have up to 820k.

   Three things made the lift fast enough, and each is only a search aid: none
   is trusted, because the final check is exact.
   - **Float pre-pass.** The working set is first grown with cheap HiGHS duals,
     so usually one exact solve suffices.
   - **Warm start.** SoPlex starts from HiGHS's optimal basis. On 2LZM χ≤3
     that took the exact solve from a 2-hour timeout to ~80 s.
   - **GMP rationals** (`gmpy2`) in the lift. Exact duals from a warm-started
     basis can be enormous: dual files of 70–310 MB.

**Result, each from running the checker on the published instance and
certificate: 42/42 exact on the reduced instances, 38/42 on the unpruned
instances, with gap exactly 0 in every case.** The four not yet lifted are
2LZM χ≤3, 1QOP χ≤3, 1GAI χ≤3 and 3PGK χ≤3, the instances with the most
triplets. Their exact certificates are on the reduced instance, so DEE is still
in their chain. The run data is in
`experiments/lp_bound.*.jsonl`, `tighten.jsonl`, `exact_certificate.R.jsonl` and
`lift_certificate.*.jsonl`. The reduced certificates
include all twelve instances CROWN could not certify: 1UBQ χ≤3, the treewidth-31 case above, and
1GAI χ≤3, where CROWN's best answer was 867.8 kcal/mol too high.

| instance | rotamers | LP gap (a) | triplets | exact, reduced | exact, unpruned (lifted) | exact rounds | final set |
|---|---|---|---|---|---|---|---|
| 1CRN χ≤1 | 104 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 2 | 35 / 104 |
| 1SHG χ≤1 | 194 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 57 / 194 |
| 1UBQ χ≤1 | 251 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 97 / 251 |
| 1QYS χ≤1 | 269 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 104 / 269 |
| 3CHY χ≤1 | 384 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 110 / 384 |
| 1AKI χ≤1 | 358 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 127 / 358 |
| 1A6M χ≤1 | 462 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 152 / 462 |
| 2LZM χ≤1 | 516 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 171 / 516 |
| 8ABP χ≤1 | 867 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 5 | 303 / 867 |
| 1ADE χ≤1 | 1284 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 5 | 441 / 1284 |
| 1GAI χ≤1 | 1257 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 439 / 1257 |
| 1QOP χ≤1 | 703 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 232 / 703 |
| 3PGK χ≤1 | 1188 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 1 | 404 / 1188 |
| 4AKE χ≤1 | 642 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 219 / 642 |
| 1CRN χ≤2 | 182 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 2 | 44 / 182 |
| 1SHG χ≤2 | 494 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 91 / 494 |
| 1UBQ χ≤2 | 614 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 161 / 614 |
| 1QYS χ≤2 | 578 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 133 / 578 |
| 3CHY χ≤2 | 975 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 207 / 975 |
| 1AKI χ≤2 | 850 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 201 / 850 |
| 1A6M χ≤2 | 1197 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 287 / 1197 |
| 2LZM χ≤2 | 1278 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 5 | 309 / 1278 |
| 8ABP χ≤2 | 2094 | **0.192** | 108 | ✔ gap 0 | ✔ gap 0 | 1 | 698 / 2094 |
| 1ADE χ≤2 | 3054 | **0.336** | 55 | ✔ gap 0 | ✔ gap 0 | 5 | 816 / 3054 |
| 1GAI χ≤2 | 2787 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 677 / 2787 |
| 1QOP χ≤2 | 1711 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 5 | 378 / 1711 |
| 3PGK χ≤2 | 2862 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 1 | 678 / 2862 |
| 4AKE χ≤2 | 1569 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 6 | 418 / 1569 |
| 1CRN χ≤3 | 254 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 46 / 254 |
| 1SHG χ≤3 | 872 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 122 / 872 |
| 1UBQ χ≤3 | 1127 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 165 / 1127 |
| 1QYS χ≤3 | 920 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 175 / 920 |
| 3CHY χ≤3 | 1614 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 239 / 1614 |
| 1AKI χ≤3 | 1399 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 3 | 236 / 1399 |
| 1A6M χ≤3 | 2142 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 4 | 504 / 2142 |
| 2LZM χ≤3 | 2214 | **0.22** | 25 | ✔ gap 0 | not run | – | – |
| 8ABP χ≤3 | 3651 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 1 | 726 / 3651 |
| 1ADE χ≤3 | 5106 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 6 | 988 / 5106 |
| 1GAI χ≤3 | 4002 | **0.00117** | 196 | ✔ gap 0 | not run | – | – |
| 1QOP χ≤3 | 2800 | **0.967** | 24 | ✔ gap 0 | not run | – | – |
| 3PGK χ≤3 | 4986 | **0.152** | 63 | ✔ gap 0 | not run | – | – |
| 4AKE χ≤3 | 2874 | ≤1e-11 | 0 | ✔ gap 0 | ✔ gap 0 | 5 | 594 / 2874 |

"exact rounds" counts exact SoPlex solves in the lift; "final set" is the
working set against all rotamers. Instances lifted before the floating-point
pre-pass existed needed more exact rounds.

**What the certificate is about, precisely.** It is about the energy arrays
**as loaded by the verifier** from the published instance files
(`packing/instance_io.py`: one `.npz` per instance with the exact float64
tables, SHA-256 in `experiments/ARTIFACTS.sha256`). The certificate binds to
them by a hash over their exact bytes. PDB parsing, rotamer building and the
force-field code are upstream of that and not certified. The instances are
published as files, rather than recomputed, because the energy function is
**not bit-reproducible across machines**. The
same PDB file gave different last bits on a second machine, and so a different
instance hash. A certificate says nothing about any other precision, force field
or machine, and nothing about whether the energy function is any good.

### What this changed about the earlier claims

- **CROWN's own verifier is tolerance-based.** `crown/verify.py` treats
  `|E − LB| ≤ max(1e-5, 1e-5·|E|)` as optimal. At these energies (1e4–8e4
  kcal/mol) that is 0.1–0.8 kcal/mol. So the 30 earlier "certified" optima were
  certified only to that tolerance, not exactly. Their exact agreement with
  `toulbar2` was the real evidence that they were right, and the exact
  certificates now confirm all 30.
- **Speed is still `toulbar2`'s.** Exact certification here costs seconds to
  minutes per instance: SoPlex ≤345 s on the reduced LPs, a median of 2.4 s.
  `toulbar2` stays under a second. The claim remains the artifact.

### Four things that went wrong on the way

These are kept on purpose; each would have produced a wrong or unverifiable
certificate.

- **SoPlex without GMP is not exact.** Built against Boost's rationals only, its
  "exact" mode prints `Cannot set optimality tolerance to small value 0 without
  GMP - using 1e-10` and still reports `optimal`. Its duals then missed the
  exact LP optimum by ~1e-14 in 11 of 12 random tests.
  - `tests/test_certify.py` demands strong duality in exact arithmetic
    (certificate LB == exact primal LP value), which caught this.
  - `scripts/build_soplex.sh` refuses to build without GMP, and the driver
    treats any SoPlex warning as fatal.
- **SoPlex's MPS reader caps lines at 256 characters.** The force field
  produces energies as small as 2.6e-268, whose exact rationals run to ~300
  digits, so the LPs are written in LP format (8190-character lines) instead.
- **The first checker was forgeable.** Two spellings of one pair key (`"0,1"` and
  `"00,1"`) inside a triplet message were added to the pair twice but subtracted
  once, so a crafted certificate could "prove" a non-optimal assignment. An
  adversarial review found it. Keys must now be canonical, shapes exact and
  duplicates are rejected, and the forgery is a regression test.
- **The instances are not reproducible across machines** (above). The
  certified objects are therefore the published instance files, not a
  recomputation.

## Correctness

Every claim above rests on checks that run in CI, because encoding bugs are the
failure mode in this area:

- **Round-trip**: for both encodings, on random instances, `QUBO(x) + const`
  equals the MRF energy of `decode(x)`, over many random assignments.
- **DEE safety**: brute-force global optimum before pruning equals brute-force
  global optimum after pruning. DEE safety is a theorem; this tests the
  implementation.
- **End-to-end lift**: each certified solution on the pruned instance is mapped
  back to *original* rotamer indices and re-evaluated on the *original*
  instance, asserting equality. This is what makes "a certificate on the pruned
  problem is a certificate on the original" a verified statement rather than an
  asserted one.
- **Agreement with an independent exact solver**: all 30 CROWN-certified optima
  equal `toulbar2`'s, to floating-point equality, on the same 42-instance set.
- **Certificate checker** (`tests/test_certify.py`):
  - Soundness against brute force: LB ≤ the true minimum for *random* messages,
    not just good ones.
  - A suboptimal assignment is never certified.
  - Tampering is rejected: a one-ulp energy change, out-of-range or non-integer
    assignments, truncated messages, messages on absent pairs, and the
    duplicate-key forgery.
  - With SoPlex, exact strong duality: LB equals the exact primal LP value, and
    the gap is exactly 0 when the LP is integral.
- **Every stored certificate re-verifies from the published files alone**
  (`experiments/verify_certificate.py`), in exact arithmetic, with no solver.

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/fetch_data.sh
.venv/bin/python tests/test_packing.py
NROT=4 .venv/bin/python experiments/encoding_comparison.py
.venv/bin/python experiments/size_scaling.py         # the size / treewidth ladder
.venv/bin/python experiments/toulbar2_benchmark.py   # needs pytoulbar2

# exact certificates
bash scripts/build_soplex.sh                          # SoPlex 8.1.0 + GMP, into third_party/
.venv/bin/python tests/test_certify.py
.venv/bin/python experiments/lp_bound.py              # rungs (a) and (b), floating point
SCOPE=red .venv/bin/python experiments/exact_certificate.py
.venv/bin/python experiments/lift_certificate.py      # reduced -> unpruned
.venv/bin/python experiments/export_instances.py      # exact instance files + hash check
```

To re-check the published certificates, you need neither SoPlex nor any other
solver:

```bash
bash scripts/fetch_artifacts.sh certs-2026-09-30     # certificates + exact instances, sha256-checked
.venv/bin/python experiments/verify_certificate.py
```

## Honest scope

**Provably-optimal side-chain packing is not new.** DEE/A\*, ILP formulations,
and weighted-CSP solvers such as `toulbar2` have solved these instances exactly
for years, and on much larger proteins than these. Nothing here beats them on
speed or size — and that is now **measured, not assumed**: `toulbar2` solves all
42 benchmark instances in under a second each. The CROWN pipeline certified 30
of them and was a median 11× slower. The exact LP certificates now cover all 42,
but they take seconds to minutes each. See
[Measured against `toulbar2`](#measured-against-toulbar2).

What is different is the **artifact**: an independently verifiable certificate.
A third party — a reviewer, a collaborator, a contract counterparty — can check
optimality from the certificate alone, by bounded arithmetic, without trusting or
rerunning the solver. With the LP-dual certificates, that check is exact
rational arithmetic over published instance files: no solver, no tolerance.

Further limits, stated plainly:

- The energy function is **crude** — counted steric clashes and contacts on heavy
  atoms, not a force field. **Biological accuracy is not claimed.**
- Rotamers come from a **staggered multi-χ grid**, not a real library. Mean
  arity reaches ~24 at χ≤3, but there are no rotamer priors, so a
  low-probability conformation costs the same as a common one.
- Proteins run to **369 flexible residues**, and all 42 instances now certify
  exactly (38 without DEE in the chain). But that rests on the LP relaxation (plus triplets) being tight on
  these instances, which is an observation, not a theorem. A more frustrated
  energy (item 7 below) could open real gaps that triplets do not close.
- A certificate is about the **stored doubles** in the published instance
  files. The energy function is not bit-reproducible across machines, so a
  recomputed instance may not be the certified one.
- A certificate says "optimal **for this energy function**." Model error is
  untouched — and in protein design model error is usually the binding
  constraint.

That last point is the actual value rather than a caveat. Optimisation error and
model error are confounded in every annealing pipeline: a bad structure could
mean a bad energy function or a stuck sampler, and you cannot tell which. **A
certificate removes one of the two unknowns**, which is what lets you attribute a
failure to the model.

## Next

1. ~~Real energy function~~ — **done.** DEE prunes 47–75% under molecular
   mechanics and certification is unchanged.
2. ~~Real rotamer libraries~~ — **done, with a caveat.** Multi-χ generation
   reaches mean arity ~24; DEE pruning *improves* to 81–90%. But 2/15 instances
   now exceed the certification ceiling at ~90–105 residual variables.
3. ~~Shrink the residual core~~ — **done.** Split DEE takes 13/15 → **15/15**,
   cutting survivors by 28–47% on the hard cases and running faster.
4. ~~Larger proteins~~ — **done, and the premise was wrong.** Certification
   reaches **369 flexible residues**, 12× the original set, and the "~90
   variable ceiling" turned out to be an artifact of a small sample: a
   251-variable core certifies while a 93-variable one does not. There is no
   size ceiling. See
   [Larger proteins](#larger-proteins-the-ceiling-is-not-a-size-ceiling).
5. ~~Attack treewidth directly~~ — **done, by not attacking it.** Pair-split
   DEE (`split=2`) was an honest negative: 1UBQ χ≤3 went from 93 to 87
   variables and 3CHY from 98 to 96, and neither crossed the line. What worked
   was replacing exact core solving with an exact **LP-dual certificate**: the
   local-polytope LP plus triplet clusters is tight on all 42 instances, and
   SoPlex plus a `Fraction` checker makes that exact. Treewidth no longer
   matters. See [Exact certificates for all 42](#exact-certificates-for-all-42-the-relaxation-not-the-core).
6. **A real library** (Dunbrack, with rotamer priors) to replace the staggered
   grid, which would also let low-probability rotamers be pruned on prior.
7. **Electrostatics and proper solvation** (EEF1 or GB) — where the landscape
   gets genuinely frustrated and persistency arguments are most likely to weaken.
8. ~~Benchmark against `toulbar2`~~ — **done, and it is a loss.** `toulbar2`
   solves 42/42 in under a second each; we certify 30/42 and are a median 11×
   slower where we succeed. The prediction that it would solve what we cannot
   was correct. The benchmark did, however, confirm all 30 certified optima
   independently. See
   [Measured against `toulbar2`](#measured-against-toulbar2).

9. **Harder instances for the relaxation.** All 42 are LP-tight after at most
   one round of triplets. Frustrated energies (electrostatics, item 7) or
   design problems with more rotamers per residue are where larger cycles or
   the SDP rung would first be needed. Finding where the LP stops being tight
   is the next honest test of this certificate.
10. **Lift the last four.** 2LZM χ≤3, 1QOP χ≤3, 1GAI χ≤3 and 3PGK χ≤3 still
    trust DEE. Their triplet clusters cover the whole working set at three
    residues, so each table grows as the cube of the rotamers added, and the
    LP slows every round. Restricting triplets to the rotamers that matter, or
    lifting triplet messages the way pair messages are lifted, should fix it.

## Licence

MIT. CROWN is Apache-2.0. PDB coordinates are from RCSB and are fetched, not
redistributed.
