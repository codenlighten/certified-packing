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
  falls back to exact core solving, which is exponential in treewidth. It holds
  to ~90 variables and fails at 100–120.

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
- **4 rotamers per residue** via χ₁ rotation; real libraries carry 10–100 and
  include χ₂₊.
- These are **small proteins** (9–30 flexible residues).
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
2. **Real rotamer libraries** — Dunbrack backbone-dependent, χ₁+χ₂, 10–100
   rotamers per residue. This is the real scaling test, and 1E0Q's 47% pruning
   rate is the warning sign to watch.
3. **Electrostatics and proper solvation** (EEF1 or GB), which is where the
   energy landscape gets genuinely frustrated.
4. **Larger proteins**, and a benchmark against `toulbar2` for time-to-optimal.
5. **Split DEE** and stronger elimination criteria for when Goldstein stalls —
   directly relevant to 1E0Q.

## Licence

MIT. CROWN is Apache-2.0. PDB coordinates are from RCSB and are fetched, not
redistributed.
