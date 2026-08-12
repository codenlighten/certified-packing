#!/usr/bin/env bash
# Fetch PDB coordinates from RCSB. Nothing is redistributed in this repository.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/pdb
for p in 1VII 1FME 2I9M 1E0Q 1L2Y; do
  echo "==> $p"
  curl -fSL --retry 3 -o "data/pdb/$p.pdb" "https://files.rcsb.org/download/$p.pdb"
done
echo "Done."
