"""Exact storage of an MRF instance: the stored doubles, bit for bit.

The energy function is not bit-reproducible across machines (numpy's SIMD
kernels differ, so the same PDB file gives doubles that differ in the last
bits), and a certificate is about the exact stored doubles. So the instances a
certificate refers to are published as files, and the verifier checks the
certificate against the file, not against a recomputation.

Format: one compressed .npz per instance with
  arity        int64[n]      rotamers per residue
  self         float64[sum]  Eself, concatenated
  pairs        int64[m, 2]   pair keys (i < j), in sorted order
  pair_data    float64[...]  Epair tables, row-major, concatenated in `pairs` order
  keep         int64[...]    Split-DEE kept rotamers per residue, concatenated
  keep_len     int64[n]
"""

import numpy as np


def save(path, Eself, Epair, keep=None):
    arity = np.array([len(e) for e in Eself], dtype=np.int64)
    pairs = sorted(Epair)
    data = {
        "arity": arity,
        "self": np.concatenate([np.asarray(e, dtype=np.float64) for e in Eself]),
        "pairs": np.array(pairs, dtype=np.int64).reshape(-1, 2),
        "pair_data": (np.concatenate([np.asarray(Epair[k], dtype=np.float64).ravel()
                                      for k in pairs]) if pairs else np.zeros(0)),
    }
    if keep is not None:
        data["keep"] = np.concatenate([np.asarray(k, dtype=np.int64) for k in keep])
        data["keep_len"] = np.array([len(k) for k in keep], dtype=np.int64)
    np.savez_compressed(path, **data)


def load(path):
    """Returns (Eself, Epair, keep or None), with exactly the stored doubles."""
    d = np.load(path)
    arity = d["arity"]
    off = np.concatenate([[0], np.cumsum(arity)])
    Eself = [d["self"][off[i]:off[i + 1]].copy() for i in range(len(arity))]
    Epair, p = {}, 0
    for i, j in d["pairs"]:
        i, j = int(i), int(j)
        size = int(arity[i] * arity[j])
        Epair[(i, j)] = d["pair_data"][p:p + size].reshape(arity[i], arity[j]).copy()
        p += size
    keep = None
    if "keep" in d:
        ko = np.concatenate([[0], np.cumsum(d["keep_len"])])
        keep = [[int(v) for v in d["keep"][ko[i]:ko[i + 1]]] for i in range(len(arity))]
    return Eself, Epair, keep
