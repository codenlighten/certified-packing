"""Lift a DEE-reduced certificate to the UNPRUNED instance, so DEE leaves the chain.

The exact certificates from `exact_certificate.py` (SCOPE=red) prove optimality
on the Split-DEE-reduced MRF. That leaves DEE, a floating-point procedure, in the
trust chain: if it ever pruned the true optimum, the reduced proof would be about
the wrong problem. Solving the unpruned LP exactly removes it, but SoPlex needs
minutes for ~40k columns, and the chi<=3 instances have up to 820k.

Lifting needs no solver. It keeps every message on the kept rotamers and gives
each PRUNED rotamer the largest message the pair terms allow, in closed form and
in exact rationals. For pair (i,j) with term minimum m (from the reduced
certificate):

  1. pruned b at j:  d_j(b) = min over KEPT a of [P(a,b) - d_i(a)] - m
  2. pruned a at i:  d_i(a) = min over ALL b  of [P(a,b) - d_j(b)] - m

where P = Epair + triplet messages. Every entry of the pair term is then >= m,
so no pair minimum moves. A triplet message entry that touches a pruned rotamer
is set to -C, with C large enough that no triplet minimum moves either. What is
left is one inequality per pruned rotamer: its node term must not fall below the
node's minimum. DEE says pruned rotamers are dominated, so this usually holds,
but nothing here assumes it: packing/certify.check() re-derives everything on the
full instance, and if lifting fails for some rotamer the certificate simply does
not verify. Each pair is lifted in whichever orientation leaves the smaller node
deficit; the orientation is a choice, not a trust assumption.

An exact dual of the REDUCED LP need not extend, even when the full LP is tight,
so a few pruned rotamers can still fail. Those are added back to the working
set, the (still small) LP is re-solved exactly by SoPlex, and lifting repeats:
row-and-column generation, with DEE's kept set only as the starting point.

The result is a certificate for the unpruned instance, checked on the unpruned
instance. The chain is then only: these messages add up, in exact arithmetic.
"""

import sys
import os
import gzip
import json
import time
from fractions import Fraction
from itertools import product

try:        # GMP rationals: exact like Fraction, much faster on large ones
    from gmpy2 import mpq
except ImportError:
    mpq = Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from packing.core import parse_residues, build_rotamers, build_instance, dee_prune
from packing import certify
from lp_bound import reduce_keep

PROTEINS = os.environ.get(
    "PROTEINS",
    "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM,8ABP,1ADE,1GAI,1QOP,3PGK,4AKE",
).split(",")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
ENERGY = os.environ.get("ENERGY", "mm")
CERT_DIR = os.environ.get("CERT_DIR", "experiments/certificates")
OUT = os.environ.get("OUT", "experiments/lift_certificate.jsonl")
WARM = os.environ.get("WARM", "1") == "1"   # HiGHS basis as SoPlex warm start
PREPASS = os.environ.get("PREPASS", "1") == "1"   # grow the working set with float duals first
F = mpq


def reduced_terms(Er, Epr, cert):
    """Exact reparameterised node/pair/triplet tables and their minima for the
    REDUCED certificate (same algebra as certify.lower_bound)."""
    n = len(Er)
    node = [[F(float(v)) for v in e] for e in Er]
    delta = {tuple(map(int, k.split(","))): ([F(v) for v in d[0]], [F(v) for v in d[1]])
             for k, d in cert["delta"].items()}
    eta = [{tuple(map(int, k.split(","))): [[F(v) for v in r] for r in rows]
            for k, rows in tm.items()} for tm in cert["eta"]]
    for (i, j), (di, dj) in delta.items():
        for a, v in enumerate(di):
            node[i][a] += v
        for b, v in enumerate(dj):
            node[j][b] += v
    pair_min = {}
    for (i, j), M in Epr.items():
        di, dj = delta.get((i, j), ([F(0)] * len(Er[i]), [F(0)] * len(Er[j])))
        extra = [tm[(i, j)] for tm in eta if (i, j) in tm]
        pair_min[(i, j)] = min(
            F(float(M[a][b])) - di[a] - dj[b] + sum((T[a][b] for T in extra), F(0))
            for a in range(len(Er[i])) for b in range(len(Er[j])))
    tri_min = []
    for (i, j, k), tm in zip(cert["triplets"], eta):
        tri_min.append(min(-(tm[(i, j)][a][b] + tm[(j, k)][b][c] + tm[(i, k)][a][c])
                           for a, b, c in product(range(len(Er[i])), range(len(Er[j])),
                                                  range(len(Er[k])))))
    return delta, eta, [min(t) for t in node], pair_min, tri_min


def lift(Eself, Epair, keep, cert_red, Er, Epr, tol=0):
    n = len(Eself)
    delta_r, eta_r, node_min, pair_min, tri_min = reduced_terms(Er, Epr, cert_red)
    pos = [{a: r for r, a in enumerate(keep[i])} for i in range(n)]   # full -> reduced
    kept = [set(keep[i]) for i in range(n)]

    # Triplet messages on the full label sets: kept entries copied, any entry
    # touching a pruned rotamer set to -C (C chosen per triplet, below).
    eta_f = []
    for t, (tm, tmin) in enumerate(zip(eta_r, tri_min)):
        emax = max((abs(v) for T in tm.values() for r in T for v in r), default=F(0))
        C = 2 * emax + abs(tmin) + 1
        full = {}
        for (p, q), T in tm.items():
            full[(p, q)] = [[(T[pos[p][a]][pos[q][b]] if a in kept[p] and b in kept[q] else -C)
                             for b in range(len(Eself[q]))] for a in range(len(Eself[p]))]
        eta_f.append(full)

    # Full pair tables P = Epair + sum of triplet messages.
    extra = {}
    for tm in eta_f:
        for pq, T in tm.items():
            extra.setdefault(pq, []).append(T)

    node_add = [[F(0)] * len(Eself[i]) for i in range(n)]
    delta_f = {}
    for (i, j), M in Epair.items():
        ki, kj = len(Eself[i]), len(Eself[j])
        di_r, dj_r = delta_r.get((i, j), ([F(0)] * len(keep[i]), [F(0)] * len(keep[j])))
        m = pair_min[(i, j)]
        ext = extra.get((i, j), [])
        P = [[F(float(M[a][b])) + sum((T[a][b] for T in ext), F(0)) for b in range(kj)]
             for a in range(ki)]

        def orient(first_j):
            di = [di_r[pos[i][a]] if a in kept[i] else None for a in range(ki)]
            dj = [dj_r[pos[j][b]] if b in kept[j] else None for b in range(kj)]
            if first_j:   # pruned b at j against kept a, then pruned a at i against all b
                for b in range(kj):
                    if dj[b] is None:
                        dj[b] = min(P[a][b] - di[a] for a in keep[i]) - m
                for a in range(ki):
                    if di[a] is None:
                        di[a] = min(P[a][b] - dj[b] for b in range(kj)) - m
            else:
                for a in range(ki):
                    if di[a] is None:
                        di[a] = min(P[a][b] - dj[b] for b in keep[j]) - m
                for b in range(kj):
                    if dj[b] is None:
                        dj[b] = min(P[a][b] - di[a] for a in range(ki)) - m
            return di, dj

        delta_f[(i, j)] = (orient(True), orient(False))

    # Pick orientations: start with first_j for all, then greedily flip a pair
    # if that lowers the total deficit of the node inequalities it touches.
    choice = {e: 0 for e in delta_f}
    theta = [[F(float(v)) for v in e] for e in Eself]

    def node_terms(i, override=None):
        t = list(theta[i])
        for e in adj[i]:
            c = override[1] if override and override[0] == e else choice[e]
            di, dj = delta_f[e][c]
            d = di if e[0] == i else dj
            for a in range(len(t)):
                t[a] += d[a]
        return t

    def deficit(i, override=None):
        t = node_terms(i, override)
        return sum((node_min[i] - v for v in t if v < node_min[i]), F(0))

    adj = {i: [] for i in range(n)}
    for e in delta_f:
        adj[e[0]].append(e)
        adj[e[1]].append(e)
    bad = {i for i in range(n) if deficit(i) > tol}
    for _ in range(3):
        if not bad:
            break
        for e in list(delta_f):
            i, j = e
            if i not in bad and j not in bad:
                continue
            now = deficit(i) + deficit(j)
            alt = deficit(i, (e, 1 - choice[e])) + deficit(j, (e, 1 - choice[e]))
            if alt < now:
                choice[e] = 1 - choice[e]
        bad = {i for i in range(n) if deficit(i) > tol}

    cert = {"version": certify.VERSION,
            "instance": certify.instance_hash(Eself, Epair),
            "assignment": [int(keep[i][cert_red["assignment"][i]]) for i in range(n)],
            "delta": {f"{i},{j}": [[str(v) for v in delta_f[(i, j)][choice[(i, j)]][0]],
                                   [str(v) for v in delta_f[(i, j)][choice[(i, j)]][1]]]
                      for (i, j) in delta_f},
            "triplets": cert_red["triplets"],
            "eta": [{f"{p},{q}": [[str(v) for v in r] for r in T] for (p, q), T in tm.items()}
                    for tm in eta_f]}
    bad_labels = []
    for i in bad:
        t = node_terms(i)
        bad_labels += [(i, a) for a, v in enumerate(t) if v < node_min[i] - tol]
    return cert, bad_labels


def float_prepass(Eself, Epair, W, x_full, triplets, tol=1e-7, max_rounds=40):
    """Grow the working set with CHEAP floating-point duals (HiGHS) before any
    exact solve. Exact duals of these LPs can be enormous rationals (hundreds of
    MB), so each exact round is slow; float duals are small dyadic rationals.
    Lifting is run on the float certificate exactly as on an exact one, and
    rotamers falling short by more than `tol` are added. This only chooses the
    working set: the certificate that is finally checked comes from an exact
    solve on it, and nothing here is trusted."""
    from lp_bound import lp_messages
    for rnd in range(max_rounds):
        t0 = time.time()
        Er, Epr = reduce_keep(Eself, Epair, W)
        _, (delta, eta), _ = lp_messages(Er, Epr, triplets)
        t_lp = time.time() - t0
        cert = {"version": certify.VERSION, "instance": certify.instance_hash(Er, Epr),
                "assignment": [W[i].index(x_full[i]) for i in range(len(W))],
                "delta": {f"{i},{j}": [[str(mpq(float(v))) for v in di],
                                       [str(mpq(float(v))) for v in dj]]
                          for (i, j), (di, dj) in delta.items()},
                "triplets": [list(t) for t in triplets],
                "eta": [{f"{p},{q}": [[str(mpq(float(v))) for v in r] for r in M]
                         for (p, q), M in tm.items()} for tm in eta]}
        t0 = time.time()
        _, bad = lift(Eself, Epair, W, cert, Er, Epr, tol=tol)
        print(f"    float round {rnd + 1}: working set {sum(len(w) for w in W)}, "
              f"short {len(bad)}, HiGHS {t_lp:.1f}s, lift {time.time() - t0:.1f}s",
              file=sys.stderr, flush=True)
        if not bad:
            return W, rnd + 1
        for i, a in bad:
            if a not in W[i]:
                W[i] = sorted(W[i] + [a])
    return W, max_rounds


def exact_reduced(Eself, Epair, W, x_full, triplets, workdir):
    """Exact reduced certificate on working set W (SoPlex), as exact_certificate.py."""
    import exact_certificate as xc
    from lp_bound import build_lp
    Er, Epr = reduce_keep(Eself, Epair, W)
    x = [W[i].index(x_full[i]) for i in range(len(W))]
    c, A, b, edges, msg_rows, tri_rows, _ = build_lp(Er, Epr, triplets)
    path = os.path.join(workdir, "lp.lp.gz")
    xc.write_lp(path, c, A, b)
    bas = xc.highs_basis(c, A, b, os.path.join(workdir, "start.bas")) if WARM else None
    y, _, t_sp, _ = xc.soplex_dual(path, A.shape[0], A.shape[1], workdir, basis=bas)
    cert = xc.to_cert(y, Er, Epr, edges, msg_rows, triplets, tri_rows, x,
                      certify.instance_hash(Er, Epr))
    return Er, Epr, cert, t_sp


def main():
    import tempfile
    print(f"LIFT reduced -> unpruned certificates (row/column generation), "
          f"energy={ENERGY}, chi {CHIS}", flush=True)
    print(f"{'protein':<8}{'chi':<5}{'labels':<8}{'dee':<6}{'final':<7}{'iters':<7}"
          f"{'soplex s':<10}{'lift s':<9}{'check s':<9}{'gap (exact)':<14}{'optimal':<8}",
          flush=True)
    for pid in PROTEINS:
        residues = parse_residues(f"data/pdb/{pid}.pdb")
        for mc in CHIS:
            flex, backbone = build_rotamers(residues, max_chi=mc)
            Eself, Epair = build_instance(flex, backbone, ENERGY)
            _, _, keep = dee_prune(Eself, Epair, split=True)
            with gzip.open(os.path.join(CERT_DIR, f"{pid}_chi{mc}_red.json.gz"), "rt") as f:
                cred = json.load(f)
            # the reduced run's assignment and triplets seed the working set
            Er, Epr = reduce_keep(Eself, Epair, keep)
            assert cred["instance"] == certify.instance_hash(Er, Epr)
            x_full = [keep[i][cred["assignment"][i]] for i in range(len(keep))]
            triplets = [tuple(t) for t in cred["triplets"]]
            W = [list(k) for k in keep]
            pre_rounds = 0
            if PREPASS:
                W, pre_rounds = float_prepass(Eself, Epair, W, x_full, triplets)
            t_sp = t_lift = 0.0
            iters = 0
            with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as wd:
                while True:
                    iters += 1
                    if iters > 1 or W != [list(k) for k in keep]:
                        Er, Epr, cred, dt = exact_reduced(Eself, Epair, W, x_full, triplets, wd)
                        t_sp += dt
                    t0 = time.time()
                    cert, bad = lift(Eself, Epair, W, cred, Er, Epr)
                    t_lift += time.time() - t0
                    print(f"    exact round {iters}: working set {sum(len(w) for w in W)}, "
                          f"short {len(bad)}, SoPlex so far {t_sp:.0f}s, lift {time.time() - t0:.0f}s",
                          file=sys.stderr, flush=True)
                    if not bad or iters >= 40:
                        break
                    for i, a in bad:          # grow the working set, keep it sorted
                        if a not in W[i]:
                            W[i] = sorted(W[i] + [a])
            t0 = time.time()
            res = certify.check(Eself, Epair, cert)
            t_ck = time.time() - t0
            name = f"{pid}_chi{mc}_full.json.gz"
            with gzip.open(os.path.join(CERT_DIR, name), "wt") as f:
                json.dump(cert, f, separators=(",", ":"))
            row = {"protein": pid, "chi": mc, "labels": sum(len(e) for e in Eself),
                   "dee_kept": sum(len(k) for k in keep), "final_kept": sum(len(w) for w in W),
                   "iters": iters, "float_rounds": pre_rounds, "soplex_s": t_sp, "lift_s": t_lift, "check_s": t_ck,
                   "energy": str(res["energy"]), "bound": str(res["bound"]),
                   "gap": str(res["gap"]), "optimal": res["optimal"], "cert": name}
            with open(OUT, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(f"{pid:<8}{mc:<5}{row['labels']:<8}{row['dee_kept']:<6}{row['final_kept']:<7}"
                  f"{iters:<7}{t_sp:<10.1f}{t_lift:<9.1f}{t_ck:<9.1f}"
                  f"{float(res['gap']):<14.3g}{str(res['optimal']):<8}", flush=True)


if __name__ == "__main__":
    main()
