"""Does DEE + certification survive realistic rotamer counts?

Step 1 showed the pipeline survives a physics-based energy at 4 rotamers per
residue. The real question is arity: published libraries carry 10-100 rotamers
per residue, and DEE's pruning power is what the whole approach rests on.

This sweeps max_chi (chi1 -> chi1+chi2 -> chi1..chi3), which drives mean arity
from ~4 to ~25, and reports whether DEE still prunes and whether certification
survives.

The warning sign to watch is 1E0Q, which pruned only 47% under molecular
mechanics at 4 rotamers while everything else pruned 65-75%.
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (parse_residues, build_rotamers, build_instance,
                          goldstein_dee, encode, evaluate)
from crown import QUBO, crown_solve, build_certificate, verify

PROTEINS = os.environ.get("PROTEINS", "2I9M,1L2Y,1E0Q,1FME,1VII").split(",")
ENERGY = os.environ.get("ENERGY", "mm")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
VAR_CAP = int(os.environ.get("VAR_CAP", 400))

print("=" * 92)
print(f"ROTAMER SCALING — energy={ENERGY}, chi depths {CHIS}")
print("=" * 92)
print(f"{'protein':<8}{'chi':<5}{'flex':<6}{'rotamers':<10}{'mean k':<8}"
      f"{'DEE kept':<10}{'pruned':<9}{'vars':<7}{'certified':<11}{'sec':<7}")

for pid in PROTEINS:
    residues = parse_residues(f"data/pdb/{pid}.pdb")
    for mc in CHIS:
        t0 = time.time()
        flex, backbone = build_rotamers(residues, max_chi=mc)
        Eself, Epair = build_instance(flex, backbone, ENERGY)
        tot = sum(len(e) for e in Eself)
        mean_k = tot / len(Eself)

        Es, Ep, keep = goldstein_dee(Eself, Epair)
        kept = sum(len(k) for k in keep)
        pruned = 100 * (tot - kept) / tot

        Q, const, decoder, lam, nv = encode(Es, Ep, "mixed")
        if nv == 0:
            cert, kind = True, "dee-complete"
        elif nv > VAR_CAP:
            cert, kind = None, f"skipped>{VAR_CAP}"
        else:
            res = crown_solve(QUBO.from_matrix(Q))
            ver = verify(QUBO.from_matrix(Q), build_certificate(
                QUBO.from_matrix(Q), res))
            choice, valid = decoder(res.assignment)
            cert = bool(res.certified_optimal) and bool(ver.certified_optimal)
            kind = res.certificate_kind
            if valid:
                e_pruned = evaluate(choice, Es, Ep)
                lifted = evaluate([keep[i][choice[i]] for i in range(len(choice))],
                                  Eself, Epair)
                assert abs(lifted - e_pruned) < 1e-6, f"LIFT MISMATCH {pid} chi{mc}"

        print(f"{pid:<8}{mc:<5}{len(flex):<6}{tot:<10}{mean_k:<8.1f}"
              f"{kept:<10}{pruned:>5.0f}%   {nv:<7}"
              f"{('YES' if cert else 'no ') + ' ' + str(kind):<22}"
              f"{time.time()-t0:>6.1f}")
