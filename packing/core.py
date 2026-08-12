"""Side-chain packing instances, Dead-End Elimination, and QUBO encodings.

An instance is a discrete pairwise MRF:
    minimise  sum_i Eself[i][r_i]  +  sum_{i<j} Epair[(i,j)][r_i, r_j]
over a choice of one rotamer r_i per flexible residue.
"""

import numpy as np
from collections import defaultdict

CONTACT_CUT = 6.0
CLASH = 3.2
NO_CHI1 = {"GLY", "ALA", "PRO"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
CHI1_WELLS = {2: [180.0], 3: [-60.0, 60.0], 4: [-60.0, 60.0, 180.0]}


# ----------------------------------------------------------------- structure
def parse_residues(path):
    res, order, want = defaultdict(dict), [], None
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[16] not in (" ", "A"):
            continue
        ch = line[21]
        want = ch if want is None else want
        if ch != want:
            continue
        if line[76:78].strip() == "H" or line[12:16].strip().startswith("H"):
            continue
        key = (line[17:20].strip(), line[22:27])
        if key not in res:
            order.append(key)
        res[key][line[12:16].strip()] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return [(k[0], k[1], res[k]) for k in order]


def _rot(points, origin, axis, deg):
    axis = axis / np.linalg.norm(axis)
    t = np.radians(deg)
    c, s = np.cos(t), np.sin(t)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return (points - origin) @ (np.eye(3) * c + s * K
                                + (1 - c) * np.outer(axis, axis)).T + origin


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def build_rotamers(residues, nrot):
    """(flex, backbone_atoms). Rotamer 0 is always the native conformation."""
    flex, fixed = [], []
    for rn, rid, atoms in residues:
        moving = [a for a in atoms if a not in BACKBONE and a != "CB"]
        if rn in NO_CHI1 or "CB" not in atoms or "N" not in atoms or not moving:
            fixed.append(np.array([atoms[a] for a in atoms]))
            continue
        CA, CB, N = atoms["CA"], atoms["CB"], atoms["N"]
        chi = _dihedral(N, CA, CB, atoms[sorted(moving)[0]])
        base = np.array([atoms[a] for a in moving])
        confs = [base] + [_rot(base, CB, CB - CA, t - chi)
                          for t in CHI1_WELLS[nrot]]
        flex.append({"name": rn, "id": rid.strip(), "confs": confs,
                     "stem": atoms["CB"]})
        fixed.append(np.array([atoms[a] for a in atoms
                               if a in BACKBONE or a == "CB"]))
    return flex, np.vstack(fixed)


# -------------------------------------------------------------------- energy
def pair_energy(A, B):
    """Bounded steric/contact score. Counting, not squared penetration."""
    d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    return float(5.0 * np.sum(d < CLASH) - 0.1 * np.sum(d < CONTACT_CUT))


def build_instance(flex, backbone):
    n = len(flex)
    Eself = [np.array([pair_energy(f["confs"][r], backbone)
                       for r in range(len(f["confs"]))]) for f in flex]
    Epair = {}
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(flex[i]["stem"] - flex[j]["stem"]) > 14.0:
                continue
            M = np.array([[pair_energy(a, b) for b in flex[j]["confs"]]
                          for a in flex[i]["confs"]])
            if np.abs(M).max() > 1e-9:
                Epair[(i, j)] = M
    return Eself, Epair


def evaluate(choice, Eself, Epair):
    e = sum(Eself[i][choice[i]] for i in range(len(Eself)))
    for (i, j), M in Epair.items():
        e += M[choice[i], choice[j]]
    return float(e)


def brute_force(Eself, Epair):
    """Exact optimum by exhaustive search. Only for small instances."""
    import itertools
    best, arg = np.inf, None
    for combo in itertools.product(*[range(len(e)) for e in Eself]):
        v = evaluate(combo, Eself, Epair)
        if v < best:
            best, arg = v, combo
    return best, arg


# ------------------------------------------------------- Dead-End Elimination
def goldstein_dee(Eself, Epair, max_rounds=50):
    """Goldstein DEE. Provably never removes the global optimum.

    Rotamer r at residue i is dead-ending if some rotamer t satisfies

        Eself[i][r] - Eself[i][t]
          + sum_j min_s ( Epair[i,j][r,s] - Epair[i,j][t,s] )  >  0

    i.e. swapping r for t improves the energy against EVERY environment. This is
    categorical persistency -- the same idea as roof duality, applied before the
    binary encoding destroys it.
    """
    n = len(Eself)
    alive = [set(range(len(e))) for e in Eself]
    nbrs = defaultdict(list)
    for (i, j), M in Epair.items():
        nbrs[i].append((j, M, False))
        nbrs[j].append((i, M, True))

    for _ in range(max_rounds):
        removed = 0
        for i in range(n):
            if len(alive[i]) <= 1:
                continue
            for r in sorted(alive[i]):
                if len(alive[i]) <= 1:
                    break
                dead = False
                for t in sorted(alive[i]):
                    if t == r:
                        continue
                    bound = Eself[i][r] - Eself[i][t]
                    for (j, M, flip) in nbrs[i]:
                        aj = sorted(alive[j])
                        # flip means M is indexed [rot_j, rot_i]
                        diff = (M[aj, r] - M[aj, t] if flip
                                else M[r, aj] - M[t, aj])
                        bound += float(np.min(diff))
                    if bound > 1e-9:
                        dead = True
                        break
                if dead:
                    alive[i].discard(r)
                    removed += 1
        if removed == 0:
            break

    keep = [sorted(a) for a in alive]
    Es = [Eself[i][keep[i]] for i in range(n)]
    Ep = {(i, j): M[np.ix_(keep[i], keep[j])] for (i, j), M in Epair.items()}
    Ep = {k: v for k, v in Ep.items() if np.abs(v).max() > 1e-9}
    return Es, Ep, keep


# ------------------------------------------------------------------ encodings
def encode(Eself, Epair, mode="onehot"):
    """Build a QUBO. Returns (Q, const, decoder, lam, nvars).

    mode='onehot' : one variable per rotamer, lambda*(sum-1)^2 per residue
    mode='mixed'  : arity-aware -- 1 rotamer fixed (no variable), 2 rotamers get
                    a single direct bit (no penalty), >=3 fall back to one-hot
    """
    n = len(Eself)
    kind, slot, nv = [], [], 0
    for i in range(n):
        k = len(Eself[i])
        if mode == "mixed" and k == 1:
            kind.append("fix"); slot.append(())
        elif mode == "mixed" and k == 2:
            kind.append("bit"); slot.append((nv,)); nv += 1
        else:
            kind.append("hot"); slot.append(tuple(range(nv, nv + k))); nv += k

    lam = 0.0
    if any(k == "hot" for k in kind):
        reach = [float(np.abs(Eself[i]).max()) for i in range(n)]
        for (i, j), M in Epair.items():
            reach[i] += float(np.abs(M).max())
            reach[j] += float(np.abs(M).max())
        lam = 1.5 * max(reach) + 1.0

    Q, const = np.zeros((nv, nv)), 0.0

    def lin(i, r):
        """Add coefficient for 'residue i takes rotamer r' -> (const|var, w)."""
        if kind[i] == "fix":
            return ("const", None)
        if kind[i] == "bit":
            return ("bit", slot[i][0])
        return ("hot", slot[i][r])

    for i in range(n):
        if kind[i] == "fix":
            const += Eself[i][0]
        elif kind[i] == "bit":
            v = slot[i][0]
            const += Eself[i][0]
            Q[v, v] += Eself[i][1] - Eself[i][0]
        else:
            for r in range(len(Eself[i])):
                v = slot[i][r]
                Q[v, v] += Eself[i][r] - lam
                for s in range(r + 1, len(Eself[i])):
                    Q[v, slot[i][s]] += 2.0 * lam

    def put(a, b, w):
        if a == b:
            Q[a, a] += w
        else:
            Q[min(a, b), max(a, b)] += w

    for (i, j), M in Epair.items():
        ki, kj = kind[i], kind[j]
        if ki == "fix" and kj == "fix":
            const += M[0, 0]
        elif ki == "fix":
            for s in range(M.shape[1]):
                t, v = lin(j, s)
                if t == "bit":
                    if s == 0:
                        const += M[0, 0]
                    else:
                        Q[v, v] += M[0, 1] - M[0, 0]
                else:
                    Q[v, v] += M[0, s]
        elif kj == "fix":
            for r in range(M.shape[0]):
                t, v = lin(i, r)
                if t == "bit":
                    if r == 0:
                        const += M[0, 0]
                    else:
                        Q[v, v] += M[1, 0] - M[0, 0]
                else:
                    Q[v, v] += M[r, 0]
        elif ki == "bit" and kj == "bit":
            a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
            u, v = slot[i][0], slot[j][0]
            const += a
            Q[u, u] += c - a
            Q[v, v] += b - a
            put(u, v, a - b - c + d)
        elif ki == "bit":
            u = slot[i][0]
            for s in range(M.shape[1]):
                v = slot[j][s]
                Q[v, v] += M[0, s]
                put(u, v, M[1, s] - M[0, s])
        elif kj == "bit":
            v = slot[j][0]
            for r in range(M.shape[0]):
                u = slot[i][r]
                Q[u, u] += M[r, 0]
                put(u, v, M[r, 1] - M[r, 0])
        else:
            for r in range(M.shape[0]):
                for s in range(M.shape[1]):
                    if M[r, s]:
                        put(slot[i][r], slot[j][s], M[r, s])

    if any(k == "hot" for k in kind):
        const += lam * sum(1 for k in kind if k == "hot")

    def decoder(x):
        x = np.asarray(x)
        ch, ok = [], True
        for i in range(n):
            if kind[i] == "fix":
                ch.append(0)
            elif kind[i] == "bit":
                ch.append(int(x[slot[i][0]]))
            else:
                on = [r for r, v in enumerate(slot[i]) if x[v] == 1]
                ok &= len(on) == 1
                ch.append(on[0] if on else 0)
        return ch, ok

    return Q, const, decoder, lam, nv


def qubo_value(Q, x):
    x = np.asarray(x, float)
    return float(x @ np.triu(Q, 1) @ x + x @ np.diag(Q))
