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
fixes it, and the reason is structural rather than incidental.

| arm | certified | mean compression | mean vars |
|---|---|---|---|
| A — one-hot (baseline) | **3/5** | 0.000 | 74 |
| B — DEE + one-hot | **5/5** | 0.656 | 25 |
| C — DEE + arity-aware | **5/5** | 0.377 | 7 |

Per instance, 4 rotamers per flexible residue:

| protein | residues | rotamers | DEE pruned | A | B | C |
|---|---|---|---|---|---|---|
| 2I9M | 9 | 36 | 69% | ✔ | ✔ | ✔ |
| 1E0Q | 16 | 64 | 67% | ✔ | ✔ | ✔ |
| 1L2Y | 13 | 52 | 60% | ✔ | ✔ | ✔ |
| 1FME | 25 | 100 | 68% | ✘ | ✔ | ✔ |
| 1VII | 30 | 120 | 68% | ✘ | ✔ | ✔ |

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

1. **Real energy function** — swap the steric counter for a proper force-field
   term (Lennard-Jones, hydrogen bonding, solvation) and re-check that DEE still
   prunes 60%+. Nothing about the method depends on the current energy.
2. **Real rotamer libraries** — Dunbrack backbone-dependent, χ₁+χ₂, 10–100
   rotamers per residue. This is the real scaling test.
3. **Larger proteins**, and a benchmark against `toulbar2` for time-to-optimal.
4. **Split DEE** and other stronger elimination criteria when Goldstein stalls.

## Licence

MIT. CROWN is Apache-2.0. PDB coordinates are from RCSB and are fetched, not
redistributed.
