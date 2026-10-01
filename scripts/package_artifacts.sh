#!/usr/bin/env bash
# Package the certificates and the exact instances they refer to for a GitHub
# release, and write the SHA-256 manifest that is committed to the repository.
#   bash scripts/package_artifacts.sh TAG
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=${1:?usage: package_artifacts.sh TAG}
mkdir -p dist
( cd experiments && sha256sum certificates/*.json.gz instances/*.npz ) > experiments/ARTIFACTS.sha256
tar czf "dist/certificates-$TAG.tar.gz" -C experiments certificates instances ARTIFACTS.sha256
sha256sum "dist/certificates-$TAG.tar.gz"
