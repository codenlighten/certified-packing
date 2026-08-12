"""Where does certification stop? A size ladder from 32 to 228 flexible residues.

Steps 1-3 established that the pipeline survives a physics-based energy and
realistic rotamer arity, and that Split DEE closes the residual core enough to
certify all 15 (protein, chi-depth) instances in the original set. But every one
of those instances certified in under 10 seconds, which means the size ceiling
was never actually probed -- the proteins were 9-30 flexible residues.

This ladder runs the same pipeline over 32-228 flexible residues to find the
real boundary. The quantity that decides everything is not protein size but
RESIDUAL VARIABLES after Split DEE: the earlier work showed the certification
ceiling sits at roughly 90 one-hot variables, and pruning percentage on a small
instance did not predict where the wall was -- total residual size did.

Reported per instance:
  flex      flexible residues
  rot       total rotamers before pruning
  gold      survivors under Goldstein DEE alone
  split     survivors under Split DEE
  vars      residual QUBO variables under the arity-aware encoding
  certified whether CROWN certified AND an independent verify() re-checked it

Timings are split into build / dee / solve so it is clear which stage is the
one that actually runs out of road.
"""

import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (parse_residues, build_rotamers, build_instance,
                          dee_prune, encode, evaluate)
from crown import QUBO, crown_solve, build_certificate, verify

LADDER = "1CRN,1SHG,1UBQ,1QYS,1A6M,3CHY,1AKI,2LZM,8ABP"
PROTEINS = os.environ.get("PROTEINS", LADDER).split(",")
ENERGY = os.environ.get("ENERGY", "mm")
CHIS = [int(c) for c in os.environ.get("CHIS", "1").split(",")]
VAR_CAP = int(os.environ.get("VAR_CAP", 400))
SPLIT = {"0": False, "1": True, "2": 2}[os.environ.get("SPLIT", "1")]
SPLIT_LABEL = {False: "OFF (Goldstein only)", True: "ON (1 witness)",
               2: "ON (witness pairs)"}[SPLIT]

print("=" * 104)
print(f"SIZE SCALING — energy={ENERGY}, chi depths {CHIS}, "
      f"Split DEE {SPLIT_LABEL}, var cap {VAR_CAP}")
print("=" * 104)
print(f"{'protein':<8}{'flex':<6}{'chi':<5}{'rot':<7}{'mean k':<8}"
      f"{'gold':<7}{'split':<7}{'pruned':<9}{'vars':<7}{'hotv':<7}{'deg':<6}"
      f"{'core':<7}{'width':<7}{'certified':<22}{'build':<8}{'dee':<8}{'solve':<8}")

rows = []
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
        tot = sum(len(e) for e in Eself)
        t_build = time.time() - t0

        t1 = time.time()
        _, _, kg = dee_prune(Eself, Epair, split=False)
        n_gold = sum(len(k) for k in kg)
        Es, Ep, keep = dee_prune(Eself, Epair, split=SPLIT)
        kept = sum(len(k) for k in keep)
        pruned = 100 * (tot - kept) / tot
        t_dee = time.time() - t1

        # Residual-core shape. Variable COUNT turned out not to predict the
        # wall, so record what the core is actually made of: residues with
        # arity >= 3 still need a one-hot penalty clique (dense, strongly
        # non-submodular), while arity-2 residues become a single free bit.
        ar = [len(e) for e in Es]
        n_hot = sum(1 for a in ar if a >= 3)
        hot_vars = sum(a for a in ar if a >= 3)
        # Mean degree of the graph the solver actually sees: edges count only
        # when BOTH endpoints survived with arity >= 2. A residue pinned to one
        # rotamer contributes a constant, not a variable, so its edges are gone.
        live = [i for i, a in enumerate(ar) if a >= 2]
        live_set = set(live)
        edges = sum(1 for (i, j) in Ep if i in live_set and j in live_set)
        deg = 2 * edges / max(1, len(live))

        t2 = time.time()
        Q, const, decoder, lam, nv = encode(Es, Ep, "mixed")
        core, width = 0, 0
        if nv == 0:
            cert, kind = True, "dee-complete"
        elif nv > VAR_CAP:
            cert, kind = None, f"skipped>{VAR_CAP}"
        else:
            q = QUBO.from_matrix(Q)
            res = crown_solve(q)
            ver = verify(q, build_certificate(q, res))
            choice, valid = decoder(res.assignment)
            cert = bool(res.certified_optimal) and bool(ver.certified_optimal)
            kind = res.certificate_kind
            # THE quantities that matter. Roof duality fixes some variables
            # outright; whatever it cannot fix is the core, which bucket
            # elimination must solve exactly at cost exponential in the core's
            # induced TREEWIDTH -- not in the number of variables. CROWN reports
            # both, so record them rather than inferring from compression.
            core = int(res.core_size)
            width = int(res.core_width)
            if valid:
                # END-TO-END: a certificate on the pruned instance is only a
                # certificate on the original if the lift back agrees exactly.
                e_pruned = evaluate(choice, Es, Ep)
                lifted = evaluate(
                    [keep[i][choice[i]] for i in range(len(choice))],
                    Eself, Epair)
                assert abs(lifted - e_pruned) < 1e-6, f"LIFT MISMATCH {pid} chi{mc}"
        t_solve = time.time() - t2

        print(f"{pid:<8}{len(flex):<6}{mc:<5}{tot:<7}{tot/len(flex):<8.1f}"
              f"{n_gold:<7}{kept:<7}{pruned:>5.0f}%   {nv:<7}{hot_vars:<7}{deg:<6.1f}"
              f"{core:<7}{width:<7}{('YES ' if cert else 'no  ') + str(kind):<22}"
              f"{t_build:<8.1f}{t_dee:<8.1f}{t_solve:<8.1f}", flush=True)
        rows.append(dict(pid=pid, flex=len(flex), chi=mc, rot=tot, vars=nv,
                         hot_vars=hot_vars, n_hot=n_hot, deg=deg, core=core,
                         width=width, cert=cert, kind=kind, pruned=pruned,
                         t=t_build + t_dee + t_solve))

done = [r for r in rows if r["cert"] is not None]
ok = [r for r in done if r["cert"]]
bad = [r for r in done if not r["cert"]]
print("\n" + "=" * 104)
print(f"certified {len(ok)}/{len(done)} attempted"
      f"   (largest certified: {max([r['flex'] for r in ok], default=0)}"
      f" flexible residues, {max([r['vars'] for r in ok], default=0)} residual vars)")


def _span(sel, key):
    v = [r[key] for r in sel]
    return f"{min(v)}-{max(v)}" if v else "n/a"


# The discriminator. If certified and failed instances OVERLAP on a quantity,
# that quantity does not explain the wall -- which is exactly what happened to
# residual variable count.
for key in ("flex", "vars", "hot_vars", "core", "width"):
    print(f"  {key:<10} certified {_span(ok, key):<12} failed {_span(bad, key)}")
print(f"  {'route':<10} certified "
      f"{ {r['kind'] for r in ok} }   failed { {r['kind'] for r in bad} }")
