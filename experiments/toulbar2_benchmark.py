"""How does this pipeline actually place against an established exact solver?

The README has always conceded that provably-optimal side-chain packing is not
new -- DEE/A*, ILP and weighted-CSP solvers have done it for years. That
concession was never measured, which made it the load-bearing unknown in the
whole repository. This measures it.

`toulbar2` (Cooper, de Givry, Schiex et al.) is a weighted-CSP solver that takes
the SAME pairwise MRF directly, with no binary encoding, and proves optimality by
soft arc consistency plus branch and bound. It is the natural incumbent.

WHAT IS COMPARED
  tb2      toulbar2 on the full, UNPRUNED instance -- the honest incumbent
  ours     Split DEE + arity-aware QUBO + CROWN, i.e. the whole pipeline

Both arms receive the identical MRF built by `packing.core.build_instance`, so
instance construction time is reported separately and excluded from both.

TWO CAVEATS, STATED UP FRONT
  1. toulbar2 works in fixed-point costs. Energies are shifted per function to be
     non-negative and discretised at `resolution` decimal digits. To keep the
     comparison honest the returned assignment is re-evaluated under the ORIGINAL
     float energy, and that is what gets compared -- not toulbar2's internal cost.
  2. These solvers do not return the same artifact. toulbar2 proves optimality
     internally; CROWN emits a certificate a third party can re-check without
     rerunning the solver. A time column cannot express that difference, so it
     should not be read as the whole story in either direction.
"""

import sys
import os
import time
import multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (parse_residues, build_rotamers, build_instance,
                          dee_prune, encode, evaluate)
from crown import QUBO, crown_solve, build_certificate, verify
import pytoulbar2 as tb2

PROTEINS = os.environ.get(
    "PROTEINS", "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM").split(",")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
ENERGY = os.environ.get("ENERGY", "mm")
VAR_CAP = int(os.environ.get("VAR_CAP", 400))
TB2_TIMEOUT = int(os.environ.get("TB2_TIMEOUT", 300))
RESOLUTION = 6


def _tb2_worker(Eself, Epair, q):
    """Build and solve inside a child process (see solve_toulbar2)."""
    m = tb2.CFN(ubinit=1e12, resolution=RESOLUTION, vac=True)
    m.Option.showSolutions = 0
    for i, e in enumerate(Eself):
        m.AddVariable(f"x{i}", [f"r{r}" for r in range(len(e))])
    for i, e in enumerate(Eself):
        m.AddFunction([f"x{i}"], (e - float(e.min())).tolist())
    for (i, j), M in Epair.items():
        m.AddFunction([f"x{i}", f"x{j}"], (M - float(M.min())).flatten().tolist())
    out = m.Solve()
    q.put(None if out is None else list(out[0]))


def solve_toulbar2(Eself, Epair):
    """Exact optimum of the MRF via toulbar2. Returns (choice, seconds, timed_out).

    Costs are shifted per function so every entry is >= 0, which toulbar2
    requires. The returned assignment is re-scored under the original float
    energy by the caller, so the shift never has to be added back.

    This build of pytoulbar2 exposes no callable time limit, so the solve runs in
    a child process that is terminated at TB2_TIMEOUT. Solver time therefore
    includes toulbar2's own model construction, which is the fair reading: it is
    work toulbar2 must do and our arm does not.
    """
    q = mp.Queue()
    p = mp.Process(target=_tb2_worker, args=(Eself, Epair, q))
    t0 = time.time()
    p.start()
    p.join(TB2_TIMEOUT)
    if p.is_alive():
        p.terminate()
        p.join()
        return None, time.time() - t0, True
    dt = time.time() - t0
    return (q.get() if not q.empty() else None), dt, False


def solve_ours(Eself, Epair):
    """Split DEE + arity-aware QUBO + CROWN. Returns (energy, certified, sec)."""
    t0 = time.time()
    Es, Ep, keep = dee_prune(Eself, Epair, split=True)
    Q, const, decoder, lam, nv = encode(Es, Ep, "mixed")
    if nv == 0:
        choice = [keep[i][0] for i in range(len(keep))]
        return evaluate(choice, Eself, Epair), True, time.time() - t0, 0
    if nv > VAR_CAP:
        return None, None, time.time() - t0, nv
    q = QUBO.from_matrix(Q)
    res = crown_solve(q)
    ver = verify(q, build_certificate(q, res))
    choice, valid = decoder(res.assignment)
    cert = bool(res.certified_optimal) and bool(ver.certified_optimal)
    e = (evaluate([keep[i][choice[i]] for i in range(len(choice))], Eself, Epair)
         if valid else None)
    return e, cert, time.time() - t0, nv


print("=" * 100)
print(f"toulbar2 BENCHMARK — energy={ENERGY}, chi depths {CHIS}")
print("=" * 100)
print(f"{'protein':<8}{'chi':<5}{'rot':<7}{'build':<8}"
      f"{'tb2 E':<14}{'tb2 s':<9}{'ours E':<14}{'ours s':<9}"
      f"{'cert':<7}{'agree':<7}")

agree_n = tb2_n = ours_n = total = 0
for pid in PROTEINS:
    path = f"data/pdb/{pid}.pdb"
    if not os.path.exists(path):
        print(f"{pid:<8}missing — run scripts/fetch_data.sh")
        continue
    residues = parse_residues(path)
    for mc in CHIS:
        t0 = time.time()
        flex, backbone = build_rotamers(residues, max_chi=mc)
        Eself, Epair = build_instance(flex, backbone, ENERGY)
        t_build = time.time() - t0
        tot = sum(len(e) for e in Eself)

        tb_choice, t_tb, tb_to = solve_toulbar2(Eself, Epair)
        # Re-score under the ORIGINAL float energy, so discretisation cannot
        # flatter either side.
        e_tb = evaluate(tb_choice, Eself, Epair) if tb_choice else None
        e_ours, cert, t_ours, nv = solve_ours(Eself, Epair)

        total += 1
        tb2_n += e_tb is not None
        ours_n += bool(cert)
        agree = (e_tb is not None and e_ours is not None
                 and abs(e_tb - e_ours) < 1e-6)
        agree_n += agree

        print(f"{pid:<8}{mc:<5}{tot:<7}{t_build:<8.1f}"
              f"{(f'{e_tb:.3f}' if e_tb is not None else 'timeout'):<14}"
              f"{t_tb:<9.2f}"
              f"{(f'{e_ours:.3f}' if e_ours is not None else f'>{VAR_CAP}v'):<14}"
              f"{t_ours:<9.2f}"
              f"{('YES' if cert else 'no'):<7}"
              f"{('=' if agree else ('DIFFER' if e_tb is not None and e_ours is not None else '-')):<7}",
              flush=True)

print("\n" + "=" * 100)
print(f"instances {total}   toulbar2 solved {tb2_n}   we certified {ours_n}   "
      f"optima agree {agree_n}/{min(tb2_n, total)}")
