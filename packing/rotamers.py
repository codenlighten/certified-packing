"""Multi-chi rotamer generation over staggered torsion wells.

NOT the Dunbrack library. The Dunbrack backbone-dependent library requires
registration and is not freely fetchable, so this generates rotamers on a
systematic grid of staggered wells -- the physical basis that published
libraries discretise and refine. It has no backbone dependence and no rotamer
priors/probabilities.

Swapping in a real library means replacing `WELLS` and adding a probability
term; nothing else in the pipeline changes.

Chi conventions are IUPAC. The moving set for a rotation about (A -> B) is every
side-chain atom strictly beyond A in greek order, which the PDB atom naming
encodes directly: CB < CG/CG1 < CD/CD1/OD1 < NE/OE1 < CZ/NZ < NH1/OH.
"""

import numpy as np

GREEK = {"A": 0, "B": 1, "G": 2, "D": 3, "E": 4, "Z": 5, "H": 6}

# per residue: the ordered chi rotation axes (first atom, second atom)
CHI_AXES = {
    "ARG": [("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ")],
    "LYS": [("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "MET": [("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "GLU": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1")],
    "GLN": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1")],
    "ASP": [("CB", "CG"), ("CG", "OD1")],
    "ASN": [("CB", "CG"), ("CG", "OD1")],
    "ILE": [("CB", "CG1"), ("CG1", "CD1")],
    "LEU": [("CB", "CG"), ("CG", "CD1")],
    "PHE": [("CB", "CG"), ("CG", "CD1")],
    "TYR": [("CB", "CG"), ("CG", "CD1")],
    "TRP": [("CB", "CG"), ("CG", "CD1")],
    "HIS": [("CB", "CG"), ("CG", "ND1")],
    "SER": [("CB", "OG")],
    "THR": [("CB", "OG1")],
    "CYS": [("CB", "SG")],
    "VAL": [("CB", "CG1")],
}

# sp3-sp3 bonds sit in staggered wells; planar/terminal groups are sampled
# on a coarser symmetric grid
WELLS_SP3 = [-60.0, 60.0, 180.0]
WELLS_PLANAR = [-90.0, 90.0]
PLANAR_LAST = {"ASP", "ASN", "PHE", "TYR", "TRP", "HIS", "GLU", "GLN"}


def greek_rank(atom_name):
    letters = [c for c in atom_name if c.isalpha()]
    if len(letters) < 2:
        return -1                      # backbone N, C, O
    return GREEK.get(letters[1], 99)


def wells_for(resname, chi_index, n_chi_full):
    """Torsion wells for chi `chi_index` (0-based) of `resname`.

    n_chi_full is the residue's FULL chi count, not a truncated one: whether a
    torsion is planar is a property of the residue, not of how many chis this
    run happens to sample. Using the truncated count put chi1 of ASP/ASN/PHE/
    TYR/TRP/HIS/GLU/GLN into planar wells, which is wrong -- chi1 is always
    sp3-sp3.
    """
    if resname in PLANAR_LAST and chi_index == n_chi_full - 1:
        return WELLS_PLANAR
    return WELLS_SP3


def _rotate(points, origin, axis, deg):
    axis = axis / np.linalg.norm(axis)
    t = np.radians(deg)
    c, s = np.cos(t), np.sin(t)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return (points - origin) @ (np.eye(3) * c + s * K
                                + (1 - c) * np.outer(axis, axis)).T + origin


def build_side_chain_rotamers(resname, atoms, max_chi=2, max_rotamers=200):
    """Enumerate rotamers for one residue.

    `atoms` maps atom name -> (xyz, element). Returns (confs, names, elems)
    where confs[0] is always the native conformation. Only atoms beyond CB move,
    so CB itself is treated as part of the fixed frame.
    """
    all_axes = CHI_AXES.get(resname, [])
    n_chi_full = len(all_axes)
    axes = all_axes[:max_chi]
    moving_names = [a for a in atoms
                    if greek_rank(a) >= 2 and a not in ("OXT",)]
    if not axes or not moving_names:
        return None

    # every axis atom must be present (NMR/X-ray models sometimes omit atoms)
    axes = [(a, b) for (a, b) in axes if a in atoms and b in atoms]
    if not axes:
        return None

    base = np.array([atoms[a][0] for a in moving_names])
    elems = [atoms[a][1] for a in moving_names]
    ranks = np.array([greek_rank(a) for a in moving_names])

    confs = [base]
    grids = [wells_for(resname, k, n_chi_full) for k in range(len(axes))]

    def recurse(conf, k):
        if len(confs) >= max_rotamers or k >= len(axes):
            return
        a, b = axes[k]
        # axis endpoints move with earlier chis, so read them from `conf`
        pa = (conf[moving_names.index(a)] if a in moving_names
              else atoms[a][0])
        pb = (conf[moving_names.index(b)] if b in moving_names
              else atoms[b][0])
        mask = ranks > greek_rank(a)
        for w in grids[k]:
            new = conf.copy()
            new[mask] = _rotate(conf[mask], pa, pb - pa, w)
            confs.append(new)
            recurse(new, k + 1)

    recurse(base, 0)
    # de-duplicate near-identical conformations
    uniq = [confs[0]]
    for c in confs[1:]:
        if all(np.abs(c - u).max() > 0.05 for u in uniq):
            uniq.append(c)
        if len(uniq) >= max_rotamers:
            break
    return uniq, moving_names, elems
