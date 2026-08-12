#!/usr/bin/env bash
# Fetch PDB coordinates from RCSB. Nothing is redistributed in this repository.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/pdb
# First five: the original small set (9-30 flexible residues).
# Remainder: the size ladder, 32-228 flexible residues (experiments/size_scaling.py).
for p in 1VII 1FME 2I9M 1E0Q 1L2Y \
         1CRN 1SHG 1UBQ 1QYS 3CHY 1AKI 1A6M 2LZM 8ABP; do
  echo "==> $p"
  curl -fSL --retry 3 -o "data/pdb/$p.pdb" "https://files.rcsb.org/download/$p.pdb"
done
echo "Done."
