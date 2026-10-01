"""Tests for the reparameterisation certificate checker (packing/certify.py).

The checker is the one piece of code a third party has to trust, so its three
properties are tested directly against brute force:
  1. SOUND: for ANY messages, even random ones, LB <= the true minimum.
  2. EXACT when it should be: with SoPlex's rational dual the certificate's LB
     equals the exact LP optimum (strong duality, which also pins the dual sign
     convention), and an integral LP certifies the optimum with gap exactly 0.
  3. BINDING: a changed energy, a wrong assignment or a malformed message is
     rejected, and a non-optimal assignment is never certified.
"""

import sys
import os
import copy
import tempfile
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))
from packing.core import brute_force
from packing import certify
from lp_bound import lp_messages

rng = np.random.default_rng(11)
fails = []


def check(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        fails.append(name)


def random_instance(n, arities, density=0.7, frustrated=False):
    Eself = [rng.normal(0, 3, size=arities[i]) for i in range(n)]
    Epair = {}
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < density:
                M = rng.normal(0, 3, size=(arities[i], arities[j]))
                if frustrated:
                    M += 6 * rng.choice([-1, 1]) * np.eye(arities[i], arities[j])
                Epair[(i, j)] = M
    return Eself, Epair


def cert_from(msgs, triplets, Eself, Epair, x, qfun=lambda v: str(Fraction(float(v)))):
    (delta, eta) = msgs
    return {"version": certify.VERSION,
            "instance": certify.instance_hash(Eself, Epair),
            "assignment": [int(v) for v in x],
            "delta": {f"{i},{j}": [[qfun(v) for v in di], [qfun(v) for v in dj]]
                      for (i, j), (di, dj) in delta.items()},
            "triplets": [list(t) for t in triplets],
            "eta": [{f"{p},{q}": [[qfun(v) for v in r] for r in M]
                     for (p, q), M in tm.items()} for tm in eta]}


def triangles(Epair, n):
    return [(i, j, k) for i in range(n) for j in range(i + 1, n)
            for k in range(j + 1, n)
            if (i, j) in Epair and (j, k) in Epair and (i, k) in Epair]


print("SOUND: LB <= true minimum for random messages")
for trial in range(8):
    n = int(rng.integers(3, 6))
    Eself, Epair = random_instance(n, [int(rng.integers(1, 4)) for _ in range(n)])
    opt, arg = brute_force(Eself, Epair)
    tri = triangles(Epair, n)
    delta = {e: (rng.normal(0, 5, len(Eself[e[0]])), rng.normal(0, 5, len(Eself[e[1]])))
             for e in Epair}
    eta = [{(i, j): rng.normal(0, 5, (len(Eself[i]), len(Eself[j]))),
            (j, k): rng.normal(0, 5, (len(Eself[j]), len(Eself[k]))),
            (i, k): rng.normal(0, 5, (len(Eself[i]), len(Eself[k])))}
           for (i, j, k) in tri]
    cert = cert_from((delta, eta), tri, Eself, Epair, arg)
    res = certify.check(Eself, Epair, cert)
    check(f"trial {trial}: random-message LB {float(res['bound']):.3f} <= opt {opt:.3f}",
          res["bound"] <= certify.exact_energy(arg, Eself, Epair))

print("\nTIGHT: LP-dual messages certify the optimum (frustrated instances need triplets)")
for trial in range(6):
    n = int(rng.integers(4, 7))
    Eself, Epair = random_instance(n, [int(rng.integers(2, 4)) for _ in range(n)],
                                   density=0.9, frustrated=True)
    opt, arg = brute_force(Eself, Epair)
    exact_opt = certify.exact_energy(arg, Eself, Epair)
    tri = triangles(Epair, n)
    _, msgs, _ = lp_messages(Eself, Epair, tri)
    res = certify.check(Eself, Epair, cert_from(msgs, tri, Eself, Epair, arg))
    gap = res["gap"]
    check(f"trial {trial}: LB <= opt, float-message gap {float(gap):.2e}",
          res["bound"] <= exact_opt)
    if abs(gap) < 1e-9:
        # optimal assignment certified; a non-optimal one must not be
        worse = None
        import itertools
        for combo in itertools.product(*[range(len(e)) for e in Eself]):
            if certify.exact_energy(combo, Eself, Epair) > exact_opt:
                worse = combo
                break
        if worse is not None:
            bad = certify.check(Eself, Epair,
                                cert_from(msgs, tri, Eself, Epair, worse))
            check(f"trial {trial}: suboptimal assignment NOT certified",
                  not bad["optimal"])

print("\nBINDING: tampering is detected")
n = 5
Eself, Epair = random_instance(n, [3] * n, density=1.0)
opt, arg = brute_force(Eself, Epair)
_, msgs, _ = lp_messages(Eself, Epair)
cert = cert_from(msgs, [], Eself, Epair, arg)
E2 = [e.copy() for e in Eself]
E2[0][0] = np.nextafter(E2[0][0], np.inf)  # one ulp
check("one-ulp energy change breaks the instance hash",
      not certify.check(E2, Epair, cert)["hash"])
bad = copy.deepcopy(cert)
bad["assignment"][0] = 99
check("out-of-range assignment rejected", not certify.check(Eself, Epair, bad)["optimal"])
bad = copy.deepcopy(cert)
k = next(iter(bad["delta"]))
bad["delta"][k][0] = bad["delta"][k][0][:-1]
try:
    certify.check(Eself, Epair, bad)
    check("truncated message rejected", False)
except ValueError:
    check("truncated message rejected", True)
bad = copy.deepcopy(cert)
bad["delta"]["0,99"] = [["0"], ["0"]]
try:
    certify.check(Eself, Epair, bad)
    check("message on a non-existent pair rejected", False)
except ValueError:
    check("message on a non-existent pair rejected", True)

# Forgery found in review: two spellings of one pair inside a triplet message
# ("0,1" and "00,1") used to be added to the pair twice but subtracted once.
Ef = [np.array([0.0, 1.0]), np.zeros(2), np.zeros(2)]
Pf = {(0, 1): np.zeros((2, 2)), (1, 2): np.zeros((2, 2)), (0, 2): np.zeros((2, 2))}
two, zero = [["2", "2"], ["2", "2"]], [["0", "0"], ["0", "0"]]
forged = {"version": certify.VERSION, "instance": certify.instance_hash(Ef, Pf),
          "assignment": [1, 0, 0], "delta": {}, "triplets": [[0, 1, 2]],
          "eta": [{"0,1": two, "00,1": two, "0,2": zero, "1,2": zero}]}
try:
    r = certify.check(Ef, Pf, forged)
    check("duplicate-spelling triplet forgery rejected", not r["optimal"])
except ValueError:
    check("duplicate-spelling triplet forgery rejected", True)
forged2 = copy.deepcopy(forged)
forged2["eta"] = [{"0,1": two, "0,2": zero, "1,2": zero}, {"0,1": two, "0,2": zero, "1,2": zero}]
try:
    certify.check(Ef, Pf, forged2)
    check("eta without a matching triplet rejected", False)
except ValueError:
    check("eta without a matching triplet rejected", True)
forged3 = copy.deepcopy(forged)
forged3["assignment"] = [1.0, 0, 0]
check("non-integer assignment rejected",
      not certify.check(Ef, Pf, {**forged3, "triplets": [], "eta": []})["optimal"])

SOPLEX = os.environ.get("SOPLEX", os.path.join(os.path.dirname(__file__), "..",
                                               "third_party/soplex/build/bin/soplex"))
if os.path.exists(SOPLEX):
    import exact_certificate as xc
    from lp_bound import build_lp
    xc.SOPLEX = SOPLEX
    print("\nEXACT: SoPlex rational duals give gap == 0 exactly")
    for trial in range(5):
        n = int(rng.integers(4, 7))
        Eself, Epair = random_instance(n, [int(rng.integers(2, 4)) for _ in range(n)],
                                       density=0.9, frustrated=trial % 2 == 1)
        if trial == 0:  # a ~1e-268 energy, as the real force field produces
            Eself[0][0] = -2.582654524940136e-268
        opt, arg = brute_force(Eself, Epair)
        tri = triangles(Epair, n)
        c, A, b, edges, msg_rows, tri_rows, noff = build_lp(Eself, Epair, tri)
        with tempfile.TemporaryDirectory() as wd:
            mps = os.path.join(wd, "lp.lp.gz")
            xc.write_lp(mps, c, A, b)
            bas = xc.highs_basis(c, A, b, os.path.join(wd, "s.bas")) if trial % 2 else None
            y, xlp, _, _ = xc.soplex_dual(mps, A.shape[0], A.shape[1], wd, basis=bas)
        cert = xc.to_cert(y, Eself, Epair, edges, msg_rows, tri, tri_rows, arg,
                          certify.instance_hash(Eself, Epair))
        res = certify.check(Eself, Epair, cert)
        primal = sum((Fraction(float(c[j])) * v for j, v in xlp.items()), Fraction(0))
        check(f"trial {trial}: certificate LB == exact LP optimum (strong duality)",
              res["bound"] == primal)
        integral = all(v in (0, 1) for v in xlp.values())
        if integral:
            check(f"trial {trial}: integral LP -> optimum certified with gap exactly 0",
                  res["optimal"] and res["gap"] == 0)
else:
    print(f"\n(SoPlex not built at {SOPLEX}; exact-dual tests skipped)")

print(f"\n{'ALL PASS' if not fails else str(len(fails)) + ' FAILURES: ' + str(fails)}")
sys.exit(1 if fails else 0)
