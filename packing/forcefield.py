"""Energy functions for side-chain packing.

Two are provided so results can be compared directly:

  steric_count        -- the coarse counting score used in the first experiments
  molecular_mechanics -- a physics-based term: softened Lennard-Jones 6-12 with
                         real per-element radii and well depths, a distance
                         hydrogen-bond term, and a burial-based solvation term

`molecular_mechanics` is a SIMPLIFIED molecular-mechanics energy, not a validated
force field. It has no electrostatics, no torsional strain, no rotamer priors,
and its solvation term is a burial count rather than EEF1/GB. It is included to
test whether the DEE + certification pipeline survives a continuous, physically
motivated energy with realistic dynamic range -- not to predict structures.
"""

import numpy as np

# van der Waals radii (A) and well depths (kcal/mol), AMBER-like
VDW_R = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
VDW_E = {"C": 0.086, "N": 0.170, "O": 0.210, "S": 0.250}
POLAR = {"N", "O"}

LIN_FRAC = 0.6       # linearise LJ below this fraction of r_min (Rosetta-style)
HB_W = 2.0           # hydrogen-bond well depth, kcal/mol
HB_R, HB_S = 2.90, 0.35
BURY_CUT = 6.0
BURY_W = {"C": -0.05, "S": -0.05, "N": 0.10, "O": 0.10}


def _params(elems):
    r = np.array([VDW_R.get(e, 1.70) for e in elems])
    e = np.array([VDW_E.get(e, 0.086) for e in elems])
    p = np.array([e in POLAR for e in elems])
    b = np.array([BURY_W.get(e, 0.0) for e in elems])
    return r, e, p, b


def steric_count(A, ea, B, eb):
    """Coarse counting score: clashes minus contacts."""
    d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    return float(5.0 * np.sum(d < 3.2) - 0.1 * np.sum(d < 6.0))


def molecular_mechanics(A, ea, B, eb):
    """Softened LJ 6-12 + hydrogen bonding + burial solvation."""
    if len(A) == 0 or len(B) == 0:
        return 0.0
    ra, ea_, pa, ba = _params(ea)
    rb, eb_, pb, bb = _params(eb)
    d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
    d = np.maximum(d, 1e-6)

    rmin = ra[:, None] + rb[None, :]
    eps = np.sqrt(ea_[:, None] * eb_[None, :])

    # --- Lennard-Jones, linearised below LIN_FRAC * rmin so clashes stay finite
    x = rmin / d
    lj = eps * (x ** 12 - 2.0 * x ** 6)
    xl = 1.0 / LIN_FRAC
    e_lin = eps * (xl ** 12 - 2.0 * xl ** 6)
    # dE/dr at the linearisation point
    slope = eps * (-12.0 * xl ** 13 + 12.0 * xl ** 7) / rmin
    r_lin = LIN_FRAC * rmin
    lj = np.where(d < r_lin, e_lin + slope * (d - r_lin), lj)
    # ignore pairs beyond any meaningful interaction
    lj = np.where(d > 8.0, 0.0, lj)

    # --- hydrogen bonding between polar pairs
    polar = pa[:, None] & pb[None, :]
    hb = -HB_W * np.exp(-((d - HB_R) ** 2) / (2 * HB_S ** 2)) * polar

    # --- burial solvation: neighbours are (un)favourable by atom polarity
    near = (d < BURY_CUT)
    solv = (ba[:, None] * near).sum() + (bb[None, :] * near).sum()

    return float(lj.sum() + hb.sum() + solv)


ENERGIES = {"steric": steric_count, "mm": molecular_mechanics}
