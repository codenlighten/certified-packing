"""Does any encoding for k>2 rotamers preserve roof-dual persistency?

BACKGROUND
----------
Side-chain packing is a discrete pairwise MRF, so it reduces to QUBO and CROWN
can in principle certify it. In practice the reduction decides everything:

  direct binary (k=2, no penalty)  -> roof duality collapses the whole problem
                                      (compression 1.000), certified at every
                                      size tested
  one-hot       (k>2, lambda term) -> roof duality fixes NOTHING
                                      (compression 0.000); certification falls
                                      back to exact core solving, which holds to
                                      ~90 variables and fails at 100-120

The one-hot penalty lambda*(sum_r x_r - 1)^2 puts a +2*lambda coupling between
every pair of rotamers within a residue -- a dense, strongly non-submodular
clique, which is exactly what roof duality cannot fix.

HYPOTHESIS
----------
Goldstein Dead-End Elimination is categorical persistency: it prunes rotamers
that provably cannot appear in the global optimum, operating on the MRF BEFORE
any binary encoding exists to destroy the structure. If DEE reduces most
residues to one or two surviving rotamers, the arity-aware ("mixed") encoding
avoids the penalty clique for those residues and persistency should survive.

ARMS (identical instances, identical energies)
  A  onehot        -- baseline
  B  dee + onehot  -- isolates the effect of pruning alone
  C  dee + mixed   -- pruning plus arity-aware encoding

DEE safety is proved in tests/test_packing.py: the global optimum is preserved
exactly. So a certificate on the pruned instance is a certificate on the
original.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (parse_residues, build_rotamers, build_instance,
                          goldstein_dee, encode, qubo_value, evaluate)
from crown import QUBO, crown_solve, build_certificate, verify

PROTEINS = ["2I9M", "1E0Q", "1L2Y", "1FME", "1VII"]
NROT = int(os.environ.get("NROT", 4))


def run(Eself, Epair, mode, label, keep=None, origS=None, origP=None):
    Q, const, decoder, lam, nv = encode(Eself, Epair, mode)
    if nv == 0:                                   # DEE solved it outright
        choice = [0] * len(Eself)
        return dict(label=label, nv=0, kind="dee-complete", certified=True,
                    compression=1.0, energy=evaluate(choice, Eself, Epair),
                    verified=True, valid=True)
    q = QUBO.from_matrix(Q)
    res = crown_solve(q)
    ver = verify(q, build_certificate(q, res))
    choice, valid = decoder(res.assignment)
    e = evaluate(choice, Eself, Epair) if valid else float("nan")
    # END-TO-END: map the pruned solution back to ORIGINAL rotamer indices and
    # re-evaluate on the ORIGINAL instance. DEE safety is a theorem; this checks
    # the implementation of it, so the certificate really covers the full problem.
    lifted = None
    if valid and keep is not None:
        orig_choice = [keep[i][choice[i]] for i in range(len(choice))]
        lifted = evaluate(orig_choice, origS, origP)
        assert abs(lifted - e) < 1e-6, (
            f"DEE LIFT MISMATCH in {label}: pruned {e} vs original {lifted}")
    # encoding self-check on the returned assignment
    if valid:
        assert abs(qubo_value(Q, res.assignment) + const - e) < 1e-6, \
            f"ENCODING BUG in {label}"
    return dict(label=label, nv=nv, kind=res.certificate_kind,
                certified=bool(res.certified_optimal),
                compression=float(res.compression_ratio), energy=e,
                verified=bool(ver.certified_optimal), valid=valid,
                lifted=lifted)


print("=" * 84)
print(f"ENCODING COMPARISON — does DEE rescue persistency?   (NROT={NROT})")
print("=" * 84)

rows = []
for pid in PROTEINS:
    residues = parse_residues(f"data/pdb/{pid}.pdb")
    flex, backbone = build_rotamers(residues, NROT)
    Eself, Epair = build_instance(flex, backbone)
    full_rot = sum(len(e) for e in Eself)

    Es, Ep, keep = goldstein_dee(Eself, Epair)
    kept = sum(len(k) for k in keep)
    arity = [len(k) for k in keep]

    print(f"\n{pid}: {len(flex)} flexible residues, {full_rot} rotamers")
    print(f"  DEE pruned {full_rot - kept}/{full_rot} rotamers "
          f"({100*(full_rot-kept)/full_rot:.0f}%) -> arities "
          f"{ {a: arity.count(a) for a in sorted(set(arity))} }")

    out = [run(Eself, Epair, "onehot", "A onehot"),
           run(Es, Ep, "onehot", "B dee+onehot", keep, Eself, Epair),
           run(Es, Ep, "mixed", "C dee+mixed", keep, Eself, Epair)]
    for r in out:
        r["pid"] = pid
        rows.append(r)
        print(f"    {r['label']:<14} vars={r['nv']:<4} "
              f"kind={str(r['kind']):<12} certified={str(r['certified']):<6} "
              f"compression={r['compression']:.3f}  E={r['energy']:.2f}"
              + (f"  lifted-to-original E={r['lifted']:.2f} OK"
                 if r.get('lifted') is not None else ""))

    cert_e = {round(r["energy"], 6) for r in out if r["valid"] and r["certified"]}
    print(f"    all CERTIFIED arms agree: {len(cert_e) == 1}"
          + ("" if out[0]["certified"] else
             f"   (uncertified baseline A: {out[0]['energy']:.2f}, "
             f"worse by {out[0]['energy'] - min(cert_e):.2f})"))

print("\n" + "=" * 84)
print("SUMMARY — certified / total, and mean compression")
print("=" * 84)
for lab in ("A onehot", "B dee+onehot", "C dee+mixed"):
    sel = [r for r in rows if r["label"] == lab]
    print(f"  {lab:<14} certified {sum(r['certified'] for r in sel)}/{len(sel)}"
          f"   mean compression {np.mean([r['compression'] for r in sel]):.3f}"
          f"   mean vars {np.mean([r['nv'] for r in sel]):.0f}")
