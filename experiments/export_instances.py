"""Write each benchmark instance's exact doubles (and Split-DEE kept set) to
experiments/instances/, and confirm that every stored certificate's instance
hash matches the exported file. That confirmation is what makes the files the
certified objects: see packing/instance_io.py for why recomputation is not
enough.
"""

import sys
import os
import gzip
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from packing.core import parse_residues, build_rotamers, build_instance, dee_prune
from packing import certify, instance_io
from lp_bound import reduce_keep

PROTEINS = os.environ.get(
    "PROTEINS",
    "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM,8ABP,1ADE,1GAI,1QOP,3PGK,4AKE",
).split(",")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
CERT_DIR = os.environ.get("CERT_DIR", "experiments/certificates")
INST_DIR = os.environ.get("INST_DIR", "experiments/instances")


def main():
    os.makedirs(INST_DIR, exist_ok=True)
    bad = 0
    for pid in PROTEINS:
        residues = parse_residues(f"data/pdb/{pid}.pdb")
        for mc in CHIS:
            flex, backbone = build_rotamers(residues, max_chi=mc)
            Eself, Epair = build_instance(flex, backbone, "mm")
            _, _, keep = dee_prune(Eself, Epair, split=True)
            path = os.path.join(INST_DIR, f"{pid}_chi{mc}.npz")
            instance_io.save(path, Eself, Epair, keep)
            E2, P2, k2 = instance_io.load(path)
            full_h = certify.instance_hash(E2, P2)
            red_h = certify.instance_hash(*reduce_keep(E2, P2, k2))
            assert full_h == certify.instance_hash(Eself, Epair)
            for f in sorted(os.listdir(CERT_DIR)):
                if not f.startswith(f"{pid}_chi{mc}_"):
                    continue
                with gzip.open(os.path.join(CERT_DIR, f), "rt") as fh:
                    h = json.load(fh)["instance"]
                want = red_h if "_red" in f else full_h
                ok = h == want
                bad += not ok
                print(f"{f:<32}{'hash matches' if ok else 'HASH MISMATCH'}", flush=True)
    print("all certificate hashes match their exported instances" if not bad
          else f"{bad} MISMATCHES")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
