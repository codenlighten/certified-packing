#!/usr/bin/env bash
# Build SoPlex 8.1.0 (exact rational LP) into third_party/, no root.
# Needs cmake, a C++ compiler, Boost headers >= 1.70 and GMP. Without GMP,
# SoPlex's "exact" mode silently keeps a 1e-10 optimality tolerance ("Cannot
# set optimality tolerance to small value 0 without GMP"), which is NOT exact:
# tests/test_certify.py catches that. If libgmp-dev is not installed, the
# distribution's signed package is fetched with `apt download` and unpacked
# locally (Debian/Ubuntu); it must match the installed libgmp10.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p third_party && cd third_party
J=${JOBS:-4}
GMP_FLAGS=()
if [ ! -f /usr/include/x86_64-linux-gnu/gmp.h ] && [ ! -f /usr/include/gmp.h ]; then
  if [ ! -f gmp/lib/libgmp.a ]; then
    rm -rf gmp-dev gmp && mkdir -p gmp/include gmp/lib
    apt download libgmp-dev >/dev/null
    dpkg-deb -x libgmp-dev_*.deb gmp-dev
    cp gmp-dev/usr/include/gmpxx.h gmp-dev/usr/include/*/gmp.h gmp/include/
    cp gmp-dev/usr/lib/*/libgmp.a gmp-dev/usr/lib/*/libgmpxx.a gmp/lib/
    rm -rf gmp-dev libgmp-dev_*.deb
  fi
  GMP_FLAGS=(-DGMP_DIR="$PWD/gmp" -DSTATIC_GMP=on)
fi
[ -d soplex ] || git clone -q https://github.com/scipopt/soplex.git
(cd soplex && git checkout -q v8.1.0)
rm -rf soplex/build
cmake -S soplex -B soplex/build -DCMAKE_BUILD_TYPE=Release -DGMP=on "${GMP_FLAGS[@]}" \
  -DBOOST=on -DMPFR=off -DPAPILO=off -DZLIB=on > /dev/null
cmake --build soplex/build -j"$J" --target soplex > /dev/null
soplex/build/bin/soplex --version 2>/dev/null | head -1
grep -q "define SOPLEX_WITH_GMP" soplex/build/soplex/config.h || { echo "GMP not linked: exact mode would not be exact" >&2; exit 1; }
