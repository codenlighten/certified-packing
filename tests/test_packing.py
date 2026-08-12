"""Correctness tests. These run before any experiment is believed.

Two properties matter and both are easy to get wrong:
  1. every encoding must round-trip: QUBO(x) + const == MRF energy of decode(x)
  2. Goldstein DEE must never remove the global optimum
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (encode, qubo_value, evaluate, brute_force,
                          goldstein_dee, dee_prune)

rng = np.random.default_rng(7)
fails = []


def check(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        fails.append(name)


def random_instance(n, arities, density=0.6):
    Eself = [rng.normal(0, 3, size=arities[i]) for i in range(n)]
    Epair = {}
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < density:
                Epair[(i, j)] = rng.normal(0, 3, size=(arities[i], arities[j]))
    return Eself, Epair


print("ROUND-TRIP: QUBO(x) + const == MRF energy")
for trial in range(6):
    n = int(rng.integers(3, 7))
    arities = [int(rng.integers(1, 5)) for _ in range(n)]
    Eself, Epair = random_instance(n, arities)
    for mode in ("onehot", "mixed"):
        Q, const, decoder, lam, nv = encode(Eself, Epair, mode)
        ok = True
        for _ in range(40):
            choice = [int(rng.integers(0, a)) for a in arities]
            # build the binary vector this choice corresponds to
            x = np.zeros(nv, int)
            cur = 0
            for i, a in enumerate(arities):
                if mode == "mixed" and a == 1:
                    continue
                if mode == "mixed" and a == 2:
                    x[cur] = choice[i]
                    cur += 1
                else:
                    x[cur + choice[i]] = 1
                    cur += a
            back, valid = decoder(x)
            ok &= valid and back == choice
            ok &= abs(qubo_value(Q, x) + const - evaluate(choice, Eself, Epair)) < 1e-8
        check(f"trial {trial} arities={arities} mode={mode}", ok)

print("\nDEE SAFETY: pruning preserves the global optimum")
for trial in range(8):
    n = int(rng.integers(3, 6))
    arities = [int(rng.integers(2, 5)) for _ in range(n)]
    Eself, Epair = random_instance(n, arities)
    before, _ = brute_force(Eself, Epair)
    Es, Ep, keep = goldstein_dee(Eself, Epair)
    after, _ = brute_force(Es, Ep)
    pruned = sum(arities) - sum(len(k) for k in keep)
    check(f"trial {trial} arities={arities} pruned={pruned} "
          f"opt {before:.4f} -> {after:.4f}", abs(before - after) < 1e-8)

print("\nSPLIT DEE SAFETY: stronger pruning still preserves the optimum")
for trial in range(10):
    n = int(rng.integers(3, 6))
    arities = [int(rng.integers(2, 5)) for _ in range(n)]
    Eself, Epair = random_instance(n, arities)
    before, _ = brute_force(Eself, Epair)
    Eg, Pg, kg = dee_prune(Eself, Epair, split=False)
    Es, Ep, keep = dee_prune(Eself, Epair, split=True)
    after, _ = brute_force(Es, Ep)
    ng = sum(len(x) for x in kg)
    ns = sum(len(x) for x in keep)
    check(f"trial {trial} arities={arities} goldstein->{ng} split->{ns} "
          f"opt {before:.4f} -> {after:.4f}",
          abs(before - after) < 1e-8 and ns <= ng)

print("\nDEE + ENCODING: pruned instance still round-trips")
for trial in range(4):
    n = int(rng.integers(4, 7))
    arities = [int(rng.integers(2, 5)) for _ in range(n)]
    Eself, Epair = random_instance(n, arities)
    Es, Ep, keep = goldstein_dee(Eself, Epair)
    ar = [len(e) for e in Es]
    Q, const, decoder, lam, nv = encode(Es, Ep, "mixed")
    ok = True
    for _ in range(30):
        choice = [int(rng.integers(0, a)) for a in ar]
        x, cur = np.zeros(nv, int), 0
        for i, a in enumerate(ar):
            if a == 1:
                continue
            if a == 2:
                x[cur] = choice[i]; cur += 1
            else:
                x[cur + choice[i]] = 1; cur += a
        ok &= abs(qubo_value(Q, x) + const - evaluate(choice, Es, Ep)) < 1e-8
    check(f"trial {trial} arities {arities} -> {ar}", ok)

print(f"\n{'ALL PASS' if not fails else str(len(fails)) + ' FAILURES: ' + str(fails)}")
sys.exit(1 if fails else 0)
