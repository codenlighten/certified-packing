"""Reparameterisation certificates: optimality proofs checked in exact arithmetic.

A certificate for a pairwise MRF

    E(x) = sum_i Eself_i(x_i) + sum_{ij} Epair_ij(x_i, x_j)

is an assignment x* plus MESSAGES:

    delta[(i,j)] = (d_i, d_j)       d_i(a) moves energy between pair (i,j) and node i
    eta[t][(p,q)]                    for each triplet t = (i,j,k), moves energy
                                     between the triplet and each of its pairs

They rewrite E, without changing it, as a sum of terms

    node    Eself_i(a) + sum_{(i,j)} d_i(a)
    pair    Epair_ij(a,b) - d_i(a) - d_j(b) + sum_{t > ij} eta_t,ij(a,b)
    triplet -eta_t,ij(a,b) - eta_t,jk(b,c) - eta_t,ik(a,c)

For every assignment the messages cancel term by term, so E(x) equals the sum of
the terms evaluated at x. Each term is at least its own minimum, so

    LB = sum of the terms' minima  <=  E(x)   for every x.

This holds for ANY messages; nothing about how they were found (an LP solver,
floating point, a bug) can make LB wrong. It can only make LB weak. If
LB >= E(x*), then x* is a global minimum.

The checker below uses only exact rationals: Python's Fraction, or GMP's mpq
with CERTIFY_ARITH=gmpy2. It needs no solver, no floating point and no
tolerance. Energies are the stored IEEE doubles, which are exact
dyadic rationals, so the claim is about exactly the stored instance.
Triplet terms add the cycle constraints of Sontag et al. (2008); with none, the
bound is the local-polytope (Schlesinger / Wainwright) LP dual.
"""

import hashlib
import json
import os
from fractions import Fraction
from itertools import product

import numpy as np

# Exact rational type. Python's Fraction by default: slow on the very large
# rationals exact LP duals can carry, but it is the standard library and easy to
# audit. CERTIFY_ARITH=gmpy2 uses GMP's mpq instead: equally exact, much faster.
# Both are exact; they must give the same verdict.
if os.environ.get("CERTIFY_ARITH") == "gmpy2":
    from gmpy2 import mpq as Q
else:
    Q = Fraction

VERSION = "reparam-cert/1"


def instance_hash(Eself, Epair):
    """sha256 over the exact bytes of every stored energy, in canonical order."""
    h = hashlib.sha256()
    h.update(VERSION.encode())
    h.update(np.array([len(Eself), len(Epair)], dtype="<i8").tobytes())
    for e in Eself:
        a = np.ascontiguousarray(e, dtype="<f8")
        h.update(np.int64(a.size).tobytes())
        h.update(a.tobytes())
    for key in sorted(Epair):
        M = np.ascontiguousarray(Epair[key], dtype="<f8")
        h.update(np.array(key, dtype="<i8").tobytes())
        h.update(np.array(M.shape, dtype="<i8").tobytes())
        h.update(M.tobytes())
    return h.hexdigest()


def _key(k, arity):
    """Parse a canonical 'i,j' key; any other spelling ('00,1', ' 0,1') is
    rejected, so two keys can never name the same pair."""
    parts = k.split(",")
    if len(parts) != arity or any(not p.isdigit() or str(int(p)) != p for p in parts):
        raise ValueError(f"non-canonical key {k!r}")
    return tuple(int(p) for p in parts)


def _table(rows, shape, what):
    """A rational matrix of exactly the given shape."""
    if len(rows) != shape[0] or any(len(r) != shape[1] for r in rows):
        raise ValueError(f"{what}: expected shape {shape}")
    return [[_q(v) for v in r] for r in rows]


def _q(x):
    """Rational from a certificate string 'p/q' or 'p', or from a float (exact)."""
    if isinstance(x, str):
        return Q(x)
    return Q(float(x))


def exact_energy(choice, Eself, Epair):
    e = sum((Q(float(Eself[i][choice[i]])) for i in range(len(Eself))),
            Q(0))
    for (i, j), M in Epair.items():
        e += Q(float(M[choice[i], choice[j]]))
    return e


def lower_bound(Eself, Epair, cert):
    """LB(messages) in exact rationals. Every node, pair and triplet table is
    enumerated in full; nothing is skipped or filtered in floating point."""
    n = len(Eself)
    delta = {}
    for k, v in cert["delta"].items():
        key = _key(k, 2)
        if key not in Epair:
            raise ValueError(f"message on pair {key}, which is not in the instance")
        if len(v) != 2:
            raise ValueError(f"message on pair {key} must be [d_i, d_j]")
        delta[key] = v
    triplets = [tuple(t) for t in cert.get("triplets", [])]
    eta = cert.get("eta", [])
    if len(triplets) != len(eta):
        raise ValueError("one eta entry is required per triplet")

    node = [[Q(float(v)) for v in Eself[i]] for i in range(n)]
    pair = {}
    for (i, j), M in Epair.items():
        P = [[Q(float(v)) for v in row] for row in M]
        ki, kj = len(Eself[i]), len(Eself[j])
        if len(P) != ki or any(len(r) != kj for r in P):
            raise ValueError(f"pair table {(i, j)} does not match its nodes")
        if (i, j) in delta:
            di = _table([delta[(i, j)][0]], (1, ki), f"d_i on {(i, j)}")[0]
            dj = _table([delta[(i, j)][1]], (1, kj), f"d_j on {(i, j)}")[0]
            for a in range(ki):
                node[i][a] += di[a]
                for b in range(kj):
                    P[a][b] -= di[a] + dj[b]
            for b in range(kj):
                node[j][b] += dj[b]
        pair[(i, j)] = P

    lb = Q(0)
    for t, tm in zip(triplets, eta):
        if len(t) != 3 or not all(isinstance(v, int) for v in t):
            raise ValueError(f"bad triplet {t!r}")
        i, j, k = t
        if not (0 <= i < j < k < n):
            raise ValueError(f"triplet {t} not sorted or out of range")
        want = {(i, j), (j, k), (i, k)}
        parsed = {}
        for key, rows in tm.items():
            pq = _key(key, 2)
            if pq in parsed or pq not in want or pq not in pair:
                raise ValueError(f"triplet {t}: bad or duplicate message key {key!r}")
            parsed[pq] = rows
        if set(parsed) != want:
            raise ValueError(f"triplet {t} must message all three pairs")
        # parse ONCE; the same tables go into the pairs and into the triplet term
        E = {pq: _table(rows, (len(Eself[pq[0]]), len(Eself[pq[1]])),
                        f"eta {t}->{pq}") for pq, rows in parsed.items()}
        for pq, T in E.items():
            for a, r in enumerate(T):
                for b, v in enumerate(r):
                    pair[pq][a][b] += v
        lb += min(-(E[(i, j)][a][b] + E[(j, k)][b][c] + E[(i, k)][a][c])
                  for a, b, c in product(range(len(Eself[i])),
                                         range(len(Eself[j])),
                                         range(len(Eself[k]))))
    for P in pair.values():
        lb += min(min(r) for r in P)
    for t in node:
        lb += min(t)
    return lb


def check(Eself, Epair, cert):
    """Verify a certificate. Returns a dict; 'optimal' is True only if the
    instance hash matches and LB >= E(x*) in exact arithmetic."""
    if cert.get("version") != VERSION:
        raise ValueError(f"unknown certificate version {cert.get('version')!r}")
    ok_hash = cert["instance"] == instance_hash(Eself, Epair)
    x = cert["assignment"]
    ok_x = (len(x) == len(Eself)
            and all(type(v) is int for v in x)
            and all(0 <= x[i] < len(Eself[i]) for i in range(len(x))))
    ub = exact_energy(x, Eself, Epair) if ok_x else None
    lb = lower_bound(Eself, Epair, cert) if ok_hash else None
    optimal = bool(ok_hash and ok_x and lb >= ub)
    return {"hash": ok_hash, "assignment": ok_x, "energy": ub, "bound": lb,
            "gap": (ub - lb) if (ub is not None and lb is not None) else None,
            "optimal": optimal}


def dump(cert, path):
    with open(path, "w") as f:
        json.dump(cert, f, separators=(",", ":"))


def load(path):
    with open(path) as f:
        return json.load(f)
