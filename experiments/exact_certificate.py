"""Exact optimality certificates, independent of treewidth and of DEE.

`lp_bound.py` showed that the local-polytope LP, plus a few triplet clusters
where needed, is tight on every benchmark instance to ~1e-12 kcal/mol. That
residual is floating-point error in the LP solver's dual, not a real gap, but a
floating-point "gap of 1e-12" is not a proof. This script removes it.

For each instance:
  1. toulbar2 supplies the assignment x* (it is only a candidate here; its
     energy is recomputed exactly and nothing else from it is trusted).
  2. The LP, with exact rational coefficients (every stored double written as
     p/q, in LP format), goes to SoPlex in exact rational mode (Gleixner, Steffy & Wolter),
     which returns an exactly optimal rational dual.
  3. The dual becomes a reparameterisation certificate (packing/certify.py),
     checked with Fractions only: optimal iff LB >= E(x*).

SCOPE=full (the default) certifies the UNPRUNED instance, so neither DEE nor
roof duality is in the chain; SCOPE=red certifies the Split-DEE-reduced one.
Triplets are chosen as in lp_bound.tighten(), on the same scope, and only when
the plain LP is not already tight in floating point.
"""

import sys
import os
import gzip
import json
import time
import subprocess
import tempfile
import re
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from packing.core import parse_residues, build_rotamers, build_instance, dee_prune
from packing import certify
from lp_bound import reduce_keep, build_lp, tighten, toulbar2_choice

PROTEINS = os.environ.get(
    "PROTEINS",
    "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM,8ABP,1ADE,1GAI,1QOP,3PGK,4AKE",
).split(",")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
ENERGY = os.environ.get("ENERGY", "mm")
SCOPE = os.environ.get("SCOPE", "full")
SOPLEX = os.environ.get("SOPLEX", "third_party/soplex/build/bin/soplex")
CERT_DIR = os.environ.get("CERT_DIR", "experiments/certificates")
OUT = os.environ.get("OUT", f"experiments/exact_certificate.{SCOPE}.jsonl")
TIME_LIMIT = int(os.environ.get("TIME_LIMIT", 7200))


def q(x):
    """Exact p/q string for a double (or 'p' when integral)."""
    f = Fraction(float(x))
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def write_lp(path, c, A, b):
    """CPLEX-LP format, every number an exact rational p/q. (MPS is not
    usable: SoPlex caps MPS lines at 256 characters, and the force field
    produces energies as small as 1e-268, whose exact rationals run to ~300
    digits. The LP reader allows 8190-character lines.) Columns have the
    default bounds [0, inf); every row is an equality."""
    A = A.tocsr()
    with gzip.open(path, "wt") as f:
        f.write("Minimize\n obj:\n")
        for j in np.flatnonzero(c):
            f.write(f" + {q(c[j])} c{j}\n")
        f.write("Subject To\n")
        for r in range(A.shape[0]):
            terms = " ".join(f"+ {q(A.data[p])} c{A.indices[p]}"
                             for p in range(A.indptr[r], A.indptr[r + 1]))
            f.write(f" r{r}: {terms} = {q(b[r])}\n")
        f.write("End\n")


def highs_basis(c, A, b, path):
    """Optimal basis from HiGHS (floating point), written as an MPS basis file
    for SoPlex's --readbas. Only a warm start: SoPlex still solves exactly, and
    nothing in the certificate depends on it. SoPlex's own floating-point simplex
    is the bottleneck on the triplet LPs (345 s vs 12 s for HiGHS on the same
    2LZM chi<=3 LP); starting from HiGHS's basis skips it."""
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    lp = highspy.HighsLp()
    m, n = A.shape
    lp.num_col_, lp.num_row_ = n, m
    lp.col_cost_ = np.asarray(c, dtype=float)
    lp.col_lower_ = np.zeros(n)
    lp.col_upper_ = np.full(n, highspy.kHighsInf)
    lp.row_lower_ = lp.row_upper_ = np.asarray(b, dtype=float)
    Acsc = A.tocsc()
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = Acsc.indptr
    lp.a_matrix_.index_ = Acsc.indices
    lp.a_matrix_.value_ = Acsc.data
    lp.a_matrix_.num_col_, lp.a_matrix_.num_row_ = n, m
    h.passModel(lp)
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return None
    basis = h.getBasis()
    B = highspy.HighsBasisStatus.kBasic
    bcols = [j for j, st in enumerate(basis.col_status) if st == B]
    nrows_nb = [r for r, st in enumerate(basis.row_status) if st != B]
    if len(bcols) != len(nrows_nb):
        return None
    with open(path, "w") as f:
        f.write("NAME basis\n")
        for j, r in zip(bcols, nrows_nb):
            f.write(f" XL c{j} r{r}\n")
        f.write("ENDATA\n")
    return path


def soplex_dual(mps, nrows, ncols, workdir, basis=None):
    """Exact rational dual y (one Fraction per row) and primal x (dict of
    nonzero columns) from SoPlex."""
    dual = os.path.join(workdir, "dual.txt")
    primal = os.path.join(workdir, "primal.txt")
    cmd = [SOPLEX, "--readmode=1", "--solvemode=2", "-f0", "-o0",
           f"-t{TIME_LIMIT}", f"-Y={dual}", f"-X={primal}", "-v3"]
    if basis:
        cmd.append(f"--readbas={basis}")
    cmd.append(mps)
    t0 = time.time()
    run = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log = run.stdout + run.stderr
    status = re.search(r"SoPlex status\s*:\s*(.*)", log)
    status = status.group(1).strip() if status else "unknown"
    if "optimal" not in status:
        raise RuntimeError(f"SoPlex: {status}\n{log[-2000:]}")
    if "Warning" in log or "without GMP" in log:
        raise RuntimeError(f"SoPlex warned; refusing its answer\n{log[-2000:]}")
    dims = re.search(r"LP has (\d+) rows (\d+) columns", log)
    if not dims or (int(dims.group(1)), int(dims.group(2))) != (nrows, ncols):
        raise RuntimeError(f"SoPlex read a different LP: {dims and dims.group(0)}")
    # kept as SoPlex's exact rational strings: parsing hundreds of thousands of
    # very large rationals into Fractions only to print them again is slow
    y = ["0"] * nrows
    with open(dual) as f:
        for line in f:
            m = re.match(r"\s*r(\d+)\s+(-?[0-9]+(?:/[0-9]+)?)\s*$", line)
            if m:
                y[int(m.group(1))] = m.group(2)
    x = {}
    with open(primal) as f:
        for line in f:
            m = re.match(r"\s*c(\d+)\s+([-+0-9/]+)\s*$", line)
            if m:
                x[int(m.group(1))] = Fraction(m.group(2))
    return y, x, dt, status


def to_cert(y, Es, Ep, edges, msg_rows, triplets, tri_rows, x, ihash):
    delta = {f"{i},{j}": [[str(y[r]) for r in ri], [str(y[r]) for r in rj]]
             for (i, j), (ri, rj) in zip(edges, msg_rows)}
    eta = [{f"{p},{qq}": [[str(y[r]) for r in row] for row in rid]
            for (p, qq), rid in rr.items()} for rr in tri_rows]
    return {"version": certify.VERSION, "instance": ihash,
            "assignment": [int(v) for v in x], "delta": delta,
            "triplets": [list(map(int, t)) for t in triplets], "eta": eta}


def main():
    os.makedirs(CERT_DIR, exist_ok=True)
    print(f"EXACT CERTIFICATES — scope={SCOPE}, energy={ENERGY}, chi {CHIS}", flush=True)
    print(f"{'protein':<8}{'chi':<5}{'labels':<8}{'cols':<10}{'trip':<6}"
          f"{'soplex s':<10}{'check s':<9}{'gap (exact)':<14}{'optimal':<8}", flush=True)
    for pid in PROTEINS:
        residues = parse_residues(f"data/pdb/{pid}.pdb")
        for mc in CHIS:
            flex, backbone = build_rotamers(residues, max_chi=mc)
            Eself, Epair = build_instance(flex, backbone, ENERGY)
            _, _, keep = dee_prune(Eself, Epair, split=True)
            Er, Epr = reduce_keep(Eself, Epair, keep)
            ch = toulbar2_choice(Er, Epr)
            x_full = [keep[i][ch[i]] for i in range(len(ch))]
            if SCOPE == "full":
                Es, Ep, x = Eself, Epair, x_full
            else:
                Es, Ep, x = Er, Epr, ch
            ub = certify.exact_energy(x, Es, Ep)
            rounds, triplets = tighten(Es, Ep, ub)
            c, A, b, edges, msg_rows, tri_rows, noff = build_lp(Es, Ep, triplets)
            with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as wd:
                mps = os.path.join(wd, "lp.lp.gz")
                write_lp(mps, c, A, b)
                y, xlp, t_sp, status = soplex_dual(mps, A.shape[0], A.shape[1], wd)
            # If the exact LP optimum is integral, it IS an optimal assignment;
            # prefer it to toulbar2's (which works at 1e-6 resolution).
            x_lp = None
            marg = [[xlp.get(j, Fraction(0)) for j in range(noff[i], noff[i + 1])]
                    for i in range(len(Es))]
            if all(sorted(m)[-1] == 1 for m in marg):
                x_lp = [m.index(Fraction(1)) for m in marg]
            tb2_agrees = x_lp is None or certify.exact_energy(x_lp, Es, Ep) == ub
            if x_lp is not None:
                x = x_lp
            cert = to_cert(y, Es, Ep, edges, msg_rows, triplets, tri_rows, x,
                           certify.instance_hash(Es, Ep))
            t0 = time.time()
            res = certify.check(Es, Ep, cert)
            t_ck = time.time() - t0
            name = f"{pid}_chi{mc}_{SCOPE}.json.gz"
            with gzip.open(os.path.join(CERT_DIR, name), "wt") as f:
                json.dump(cert, f, separators=(",", ":"))
            row = {"protein": pid, "chi": mc, "scope": SCOPE,
                   "labels": sum(len(e) for e in Es), "cols": A.shape[1],
                   "rows": A.shape[0], "triplets": len(triplets), "lp_integral": x_lp is not None,
                   "tb2_agrees": tb2_agrees,
                   "soplex_s": t_sp, "check_s": t_ck, "status": status,
                   "energy": str(res["energy"]), "bound": str(res["bound"]),
                   "gap": str(res["gap"]), "optimal": res["optimal"],
                   "cert": name}
            with open(OUT, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(f"{pid:<8}{mc:<5}{row['labels']:<8}{row['cols']:<10}"
                  f"{len(triplets):<6}{t_sp:<10.1f}{t_ck:<9.1f}"
                  f"{float(res['gap']):<14.3g}{str(res['optimal']):<8}", flush=True)


if __name__ == "__main__":
    main()
