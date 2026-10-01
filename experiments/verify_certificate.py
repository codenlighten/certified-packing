"""Re-check stored optimality certificates without any solver.

    .venv/bin/python experiments/verify_certificate.py              # all of them
    .venv/bin/python experiments/verify_certificate.py 1UBQ_chi3_full

For each certificate this loads the exact instance it refers to from
experiments/instances/ (the stored doubles, bit for bit; see
packing/instance_io.py), checks the certificate's instance hash against it, and
runs packing/certify.check: exact rational arithmetic only, no LP
solver, no toulbar2, no CROWN, no floating-point tolerance. "optimal" means
LB >= E(x*) exactly, i.e. no assignment of the stored instance has lower energy.

Set REBUILD=1 to recompute the instances from the PDB files instead. That only
matches on hardware that reproduces the reference doubles bit for bit; the
energies are NOT bit-reproducible across machines.

What this does NOT re-check is that the instance is the one you meant: the
certificate is about the stored doubles, not about any other force field,
precision or machine.
"""

import sys
import os
import gzip
import json
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from packing.core import parse_residues, build_rotamers, build_instance, dee_prune
from packing import certify, instance_io
from lp_bound import reduce_keep

CERT_DIR = os.environ.get("CERT_DIR", "experiments/certificates")
ENERGY = os.environ.get("ENERGY", "mm")
INST_DIR = os.environ.get("INST_DIR", "experiments/instances")
REBUILD = os.environ.get("REBUILD") == "1"


def instance_for(pid, chi, scope):
    if not REBUILD:
        Eself, Epair, keep = instance_io.load(os.path.join(INST_DIR, f"{pid}_chi{chi}.npz"))
        return reduce_keep(Eself, Epair, keep) if scope == "red" else (Eself, Epair)
    flex, backbone = build_rotamers(parse_residues(f"data/pdb/{pid}.pdb"), max_chi=chi)
    Eself, Epair = build_instance(flex, backbone, ENERGY)
    if scope == "red":
        _, _, keep = dee_prune(Eself, Epair, split=True)
        return reduce_keep(Eself, Epair, keep)
    return Eself, Epair


def main(names):
    if not names:
        names = sorted(f[:-len(".json.gz")] for f in os.listdir(CERT_DIR)
                       if f.endswith(".json.gz"))
    bad = 0
    print(f"{'certificate':<26}{'energy (kcal/mol)':<22}{'gap (exact)':<14}{'optimal':<9}{'s':<6}")
    for name in names:
        pid, chi, scope = name.split("_", 2)
        scope = "red" if scope == "red" else "full"
        with gzip.open(os.path.join(CERT_DIR, name + ".json.gz"), "rt") as f:
            cert = json.load(f)
        t0 = time.time()
        Eself, Epair = instance_for(pid, int(chi[3:]), scope)
        res = certify.check(Eself, Epair, cert)
        ok = res["optimal"]
        bad += not ok
        print(f"{name:<26}{float(res['energy']) if res['energy'] is not None else float('nan'):<22.6f}"
              f"{str(res['gap']):<14}{str(ok):<9}{time.time() - t0:<6.1f}", flush=True)
    print(f"\n{len(names) - bad}/{len(names)} certificates verify")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
