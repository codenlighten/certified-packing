"""Can a relaxation bound certify what exact core solving cannot?

Twelve of the 42 benchmark instances fail to certify, and the README traces
every clean failure to core treewidth > 27: CROWN's second stage solves the
roof-duality core exactly, which is exponential in width. A relaxation bound has
no such dependence. If a lower bound LB meets the energy UB of a known
assignment, that assignment is optimal, whatever the treewidth.

This measures the cheapest rung of the ladder: the LOCAL MARGINAL POLYTOPE LP
of the categorical MRF (one rotamer per residue as a hard constraint, not a
penalty). It is what RLT gives on the one-hot formulation. Its dual is a
REPARAMETERISATION: messages delta_{e->i}(a) that move energy between each pair
table and its endpoint tables. For ANY messages, whatever their source,

    LB(delta) = sum_i  min_a   [ Eself_i(a) + sum_{e ni i} delta_{e->i}(a) ]
              + sum_ij min_a,b [ Epair_ij(a,b) - delta_{ij->i}(a) - delta_{ij->j}(b) ]

is a lower bound on every assignment's energy (each assignment's energy splits
into exactly these terms, and each term is at least its minimum). So the
certificate is just the list of messages. The checker needs no LP solver, no
feasibility tolerance and no eigenvalue bound. Messages come from HiGHS's dual
values in floating point; the bound is then evaluated in EXACT rational
arithmetic over the stored doubles (every double is a dyadic rational). Solver
error can make the bound weaker, but it cannot make it wrong.

UB is the energy of toulbar2's assignment, recomputed exactly. Where LB >= UB,
the assignment is optimal for the stored instance. Otherwise the table reports
the exact gap UB - LB.

Two scopes are reported:
  red   the Split-DEE-reduced MRF. Correctness then also rests on DEE.
  full  the unpruned MRF. The chain is then only "the messages add up".
"""

import sys
import os
import time
import json
import multiprocessing as mp
from fractions import Fraction

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from packing.core import (parse_residues, build_rotamers, build_instance,
                          dee_prune, evaluate)

PROTEINS = os.environ.get(
    "PROTEINS",
    "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM,8ABP,1ADE,1GAI,1QOP,3PGK,4AKE",
).split(",")
CHIS = [int(c) for c in os.environ.get("CHIS", "1,2,3").split(",")]
ENERGY = os.environ.get("ENERGY", "mm")
FULL_CAP = int(os.environ.get("FULL_CAP", 3_000_000))  # LP columns for "full"
TIGHTEN = os.environ.get("TIGHTEN") == "1"
OUT = os.environ.get("OUT", "experiments/lp_bound.jsonl")
EPS = np.finfo(float).eps


def reduce_keep(Eself, Epair, keep):
    """The DEE-reduced MRF, keeping EVERY pair table (dee_prune drops tables
    whose entries are all below 1e-9; dropping them would change the energy)."""
    Es = [Eself[i][keep[i]] for i in range(len(Eself))]
    Ep = {(i, j): M[np.ix_(keep[i], keep[j])] for (i, j), M in Epair.items()}
    return Es, Ep


def build_lp(Es, Ep, triplets=()):
    """The local-polytope LP (plus triplet clusters) as (c, A, b, index maps).

    Columns: mu_i(a) for every node label, mu_ij(a,b) for every pair entry.
    Rows:    sum_a mu_i(a) = 1                      (dual y_i)
             sum_b mu_ij(a,b) - mu_i(a) = 0         (dual delta_{ij->i}(a))
             sum_a mu_ij(a,b) - mu_j(b) = 0         (dual delta_{ij->j}(b))
    Reduced costs are Eself_i(a) - y_i + sum delta and
    Epair_ij(a,b) - delta_i(a) - delta_j(b): exactly the reparameterisation.

    Each triplet (i,j,k) (all three pairs must be in Ep) adds columns
    mu_ijk(a,b,c) and rows sum_c mu_ijk(a,b,c) - mu_ij(a,b) = 0 for its three
    pairs (dual eta_{ijk->ij}(a,b)): the pair term gains +eta, and a triplet
    term -sum eta appears. This is Sontag et al. (2008) cycle tightening.
    """
    n = len(Es)
    noff = np.cumsum([0] + [len(e) for e in Es])
    edges = list(Ep)
    eoff = [noff[-1]]
    for (i, j) in edges:
        eoff.append(eoff[-1] + len(Es[i]) * len(Es[j]))
    ncol = eoff[-1]
    c = np.concatenate([np.concatenate(Es)] +
                       [Ep[e].ravel() for e in edges]) if edges else np.concatenate(Es)

    rows, cols, vals = [], [], []
    b = []
    r = 0
    for i in range(n):  # normalisation
        k = len(Es[i])
        rows += [r] * k
        cols += list(range(noff[i], noff[i] + k))
        vals += [1.0] * k
        b.append(1.0)
        r += 1
    msg_rows = []
    for t, (i, j) in enumerate(edges):
        ki, kj = len(Es[i]), len(Es[j])
        base = eoff[t]
        grid = base + np.arange(ki * kj).reshape(ki, kj)
        ri = []
        for a in range(ki):
            rows += [r] * (kj + 1)
            cols += list(grid[a]) + [noff[i] + a]
            vals += [1.0] * kj + [-1.0]
            b.append(0.0)
            ri.append(r)
            r += 1
        rj = []
        for bb in range(kj):
            rows += [r] * (ki + 1)
            cols += list(grid[:, bb]) + [noff[j] + bb]
            vals += [1.0] * ki + [-1.0]
            b.append(0.0)
            rj.append(r)
            r += 1
        msg_rows.append((ri, rj))
    epos = {e: t for t, e in enumerate(edges)}
    ncol = eoff[-1]
    c = list(c)
    tri_rows = []
    for (i, j, k) in triplets:
        ki, kj, kk = len(Es[i]), len(Es[j]), len(Es[k])
        tgrid = ncol + np.arange(ki * kj * kk).reshape(ki, kj, kk)
        ncol += ki * kj * kk
        c += [0.0] * (ki * kj * kk)
        rr = {}
        for pair, axis in (((i, j), 2), ((j, k), 0), ((i, k), 1)):
            t = epos[pair]
            egrid = eoff[t] + np.arange(len(Es[pair[0]]) * len(Es[pair[1]])
                                        ).reshape(len(Es[pair[0]]), len(Es[pair[1]]))
            summed = np.moveaxis(tgrid, axis, -1)  # (pair0, pair1, summed-out)
            rid = np.empty(egrid.shape, dtype=int)
            for a in range(egrid.shape[0]):
                for bb in range(egrid.shape[1]):
                    rows += [r] * (summed.shape[2] + 1)
                    cols += list(summed[a, bb]) + [egrid[a, bb]]
                    vals += [1.0] * summed.shape[2] + [-1.0]
                    b.append(0.0)
                    rid[a, bb] = r
                    r += 1
            rr[pair] = rid
        tri_rows.append(rr)
    c = np.array(c)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, ncol))
    return c, A, np.array(b), edges, msg_rows, tri_rows, noff


def lp_messages(Es, Ep, triplets=()):
    """Solve the local-polytope LP; return (LP value, messages, node marginals)."""
    n = len(Es)
    c, A, b, edges, msg_rows, tri_rows, noff = build_lp(Es, Ep, triplets)
    res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if res.status != 0:
        raise RuntimeError(res.message)
    y = res.eqlin.marginals
    msgs = {e: (y[ri], y[rj]) for e, (ri, rj) in zip(edges, msg_rows)}
    tmsgs = [{pair: y[rid] for pair, rid in rr.items()} for rr in tri_rows]
    marg = [res.x[noff[i]:noff[i + 1]] for i in range(n)]
    return float(res.fun), (msgs, tmsgs), marg


def exact_min(vals_float, terms):
    """Exact minimum of sum(terms) (a list of equally shaped float arrays),
    evaluated in rationals. Only entries whose float value is within a
    rigorous rounding bound of the float minimum can be the exact minimum,
    so only those are evaluated exactly."""
    absum = sum(np.abs(t) for t in terms)
    tol = 4 * len(terms) * EPS * float(absum.max()) + 1e-300
    cand = np.flatnonzero(vals_float.ravel() <= vals_float.min() + 2 * tol)
    flat = [t.ravel() for t in terms]
    return min(sum((Fraction(float(f[k])) for f in flat), Fraction(0))
               for k in cand)


def exact_lb(Es, Ep, messages, triplets=()):
    """LB(delta, eta) in exact rational arithmetic (see module docstring)."""
    msgs, tmsgs = messages
    node_terms = [[e] for e in Es]
    pair_extra = {e: [] for e in Ep}
    lb = Fraction(0)
    for (i, j, k), tm in zip(triplets, tmsgs):
        for pair, eta in tm.items():
            pair_extra[pair].append(eta)
        shape = (len(Es[i]), len(Es[j]), len(Es[k]))
        terms = [np.broadcast_to(-tm[(i, j)][:, :, None], shape),
                 np.broadcast_to(-tm[(j, k)][None, :, :], shape),
                 np.broadcast_to(-tm[(i, k)][:, None, :], shape)]
        lb += exact_min(terms[0] + terms[1] + terms[2], terms)
    for (i, j), M in Ep.items():
        di, dj = msgs[(i, j)]
        node_terms[i].append(di)
        node_terms[j].append(dj)
        terms = [M, np.broadcast_to(-di[:, None], M.shape),
                 np.broadcast_to(-dj[None, :], M.shape)] + pair_extra[(i, j)]
        lb += exact_min(sum(terms), terms)
    for terms in node_terms:
        lb += exact_min(sum(terms), terms)
    return lb


def exact_energy(choice, Es, Ep):
    e = sum((Fraction(float(Es[i][choice[i]])) for i in range(len(Es))), Fraction(0))
    for (i, j), M in Ep.items():
        e += Fraction(float(M[choice[i], choice[j]]))
    return e


def _tb2_worker(Es, Ep, q):
    import pytoulbar2 as tb2
    m = tb2.CFN(ubinit=1e12, resolution=6, vac=True)
    m.Option.showSolutions = 0
    for i, e in enumerate(Es):
        m.AddVariable(f"x{i}", [f"r{r}" for r in range(len(e))])
    for i, e in enumerate(Es):
        m.AddFunction([f"x{i}"], (e - float(e.min())).tolist())
    for (i, j), M in Ep.items():
        m.AddFunction([f"x{i}", f"x{j}"], (M - float(M.min())).flatten().tolist())
    out = m.Solve()
    q.put(None if out is None else list(out[0]))


def toulbar2_choice(Es, Ep, timeout=300):
    q = mp.Queue()
    p = mp.Process(target=_tb2_worker, args=(Es, Ep, q))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        p.join()
        return None
    return q.get() if not q.empty() else None


def bound(Es, Ep, ub):
    t0 = time.time()
    lpv, msgs, _ = lp_messages(Es, Ep)
    t_lp = time.time() - t0
    t0 = time.time()
    lb = exact_lb(Es, Ep, msgs)
    return {"lp": lpv, "lb": float(lb), "gap": float(ub - lb),
            "cert": lb >= ub, "t_lp": t_lp, "t_check": time.time() - t0}


def tighten(Es, Ep, ub, rounds=6, frac_tol=1e-6):
    """Rung (b): add triplet clusters on triangles touching fractional nodes
    until the exact bound meets UB or nothing is fractional."""
    trip, out = [], []
    nbr = {i: set() for i in range(len(Es))}
    for (i, j) in Ep:
        nbr[i].add(j)
        nbr[j].add(i)
    for rnd in range(rounds + 1):
        t0 = time.time()
        lpv, msgs, marg = lp_messages(Es, Ep, trip)
        t_lp = time.time() - t0
        lb = exact_lb(Es, Ep, msgs, trip)
        frac = {i for i, m in enumerate(marg) if m.max() < 1 - frac_tol}
        out.append({"round": rnd, "triplets": len(trip), "frac": len(frac),
                    "lp": lpv, "gap": float(ub - lb), "cert": lb >= ub,
                    "t_lp": t_lp})
        print(f"    round {rnd}: triplets {len(trip):<6} fractional {len(frac):<5}"
              f"gap {float(ub - lb):.3g}   lp {t_lp:.1f}s", flush=True)
        if lb >= ub or not frac:
            break
        if rnd == rounds:
            break
        have = set(trip)
        for i in sorted(frac):
            for j in nbr[i]:
                for k in nbr[i] & nbr[j]:
                    t = tuple(sorted((i, j, k)))
                    if t not in have and len(frac & set(t)) >= 2:
                        have.add(t)
                        trip.append(t)
    return out, trip


def yn(b):
    return "-" if b is None else ("Y" if b["cert"] else "n")


if __name__ == "__main__":
    print(f"LOCAL-POLYTOPE LP BOUND — energy={ENERGY}, chi {CHIS}")
    hdr = (f"{'protein':<8}{'chi':<5}{'rot':<7}{'red':<6}{'UB':<14}"
           f"{'gap red':<12}{'gap full':<12}{'cert red/full':<15}{'t_lp':<8}")
    print(hdr, flush=True)
    for pid in PROTEINS:
        path = f"data/pdb/{pid}.pdb"
        if not os.path.exists(path):
            print(f"{pid:<8}missing")
            continue
        residues = parse_residues(path)
        for mc in CHIS:
            flex, backbone = build_rotamers(residues, max_chi=mc)
            Eself, Epair = build_instance(flex, backbone, ENERGY)
            _, _, keep = dee_prune(Eself, Epair, split=True)
            Es, Ep = reduce_keep(Eself, Epair, keep)
            ch = toulbar2_choice(Es, Ep)
            full_choice = [keep[i][ch[i]] for i in range(len(ch))]
            ub = exact_energy(full_choice, Eself, Epair)
            assert ub == exact_energy(ch, Es, Ep)
            rot = sum(len(e) for e in Eself)
            red = sum(len(e) for e in Es)
            row = {"protein": pid, "chi": mc, "rot": rot, "red": red,
                   "ub": float(ub), "red_bound": bound(Es, Ep, ub)}
            ncol_full = rot + sum(M.size for M in Epair.values())
            row["full_cols"] = ncol_full
            row["full_bound"] = (bound(Eself, Epair, ub)
                                 if ncol_full <= FULL_CAP and not TIGHTEN else None)
            if TIGHTEN:
                row["tighten"], row["triplets"] = tighten(Es, Ep, ub)
            with open(OUT, "a") as f:
                f.write(json.dumps(row) + "\n")
            rb, fb = row["red_bound"], row["full_bound"]
            print(f"{pid:<8}{mc:<5}{rot:<7}{red:<6}{float(ub):<14.4f}"
                  f"{rb['gap']:<12.3g}"
                  f"{(format(fb['gap'], '.3g') if fb else 'skip'):<12}"
                  f"{yn(rb) + '/' + yn(fb):<15}"
                  f"{rb['t_lp'] + (fb['t_lp'] if fb else 0):<8.1f}", flush=True)
