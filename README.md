# certified-packing

Protein side-chain packing with a **machine-checkable proof of optimality**.

Side-chain packing — choose one rotamer per residue on a fixed backbone to
minimise a pairwise energy — is a discrete pairwise MRF. This repository poses it
as a QUBO, solves it with [CROWN](https://github.com/codenlighten/crown), and
emits a certificate a third party can verify **by bounded arithmetic, without
rerunning the solver**.

The resulting statement is categorically different from what annealing gives:

> **"This packing is the proven global minimum of the stated energy function"**
> — not "this is the best my sampler found."

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

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/fetch_data.sh
.venv/bin/python tests/test_packing.py
NROT=4 .venv/bin/python experiments/encoding_comparison.py
```

## Honest scope

**Provably-optimal side-chain packing is not new.** DEE/A\*, ILP formulations,
and weighted-CSP solvers such as `toulbar2` have solved these instances exactly
for years, and on much larger proteins than these. Nothing here beats them on
speed or size.

What is different is the **artifact**: an independently verifiable certificate.
A third party — a reviewer, a collaborator, a contract counterparty — can check
optimality from the certificate alone, by bounded arithmetic, without trusting or
rerunning the solver. That is CROWN's contribution, and this repository shows it
carries over to a real structural-biology problem.

Further limits, stated plainly:

- The energy function is **crude** — counted steric clashes and contacts on heavy
  atoms, not a force field. **Biological accuracy is not claimed.**
- Rotamers come from a **staggered multi-χ grid**, not a real library. Mean
  arity reaches ~24 at χ≤3, but there are no rotamer priors, so a
  low-probability conformation costs the same as a common one.
- Proteins run to **369 flexible residues**, but the hardest setting is not
  solved: at χ≤3 only **4/13** instances certify. Size is not the limit —
  treewidth of the residual core is.
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
5. **Attack treewidth directly** — the real binding constraint. Pair-split DEE
   (`split=2`) was tried and is an honest negative: it cut 1UBQ χ≤3 from 93 to
   87 variables and 3CHY from 98 to 96, and neither crossed the line. Pruning
   *more rotamers* is not the lever; reducing the *interaction graph* is.
   Candidates: residue-pair clustering, or a tree decomposition that solves
   high-width regions separately.
6. **A real library** (Dunbrack, with rotamer priors) to replace the staggered
   grid, which would also let low-probability rotamers be pruned on prior.
7. **Electrostatics and proper solvation** (EEF1 or GB) — where the landscape
   gets genuinely frustrated and persistency arguments are most likely to weaken.
8. **Benchmark against `toulbar2`** for time-to-optimal, to place this honestly
   against the established exact solvers. Now more pointed than before: at χ≤3
   `toulbar2` would likely solve the instances this pipeline cannot certify.

## Licence

MIT. CROWN is Apache-2.0. PDB coordinates are from RCSB and are fetched, not
redistributed.
