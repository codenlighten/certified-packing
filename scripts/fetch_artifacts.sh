#!/usr/bin/env bash
# Download the certificates and exact instances from the GitHub release and
# check every file against the committed manifest experiments/ARTIFACTS.sha256.
#   bash scripts/fetch_artifacts.sh TAG
# then: .venv/bin/python experiments/verify_certificate.py
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=${1:?usage: fetch_artifacts.sh TAG}
URL="https://github.com/codenlighten/certified-packing/releases/download/$TAG/certificates-$TAG.tar.gz"
tmp=$(mktemp -d)
curl -sSL -o "$tmp/a.tar.gz" "$URL"
tar xzf "$tmp/a.tar.gz" -C "$tmp"
# the manifest in the repository, not the one in the tarball, is authoritative
( cd "$tmp" && sha256sum -c --quiet "$OLDPWD/experiments/ARTIFACTS.sha256" )
mkdir -p experiments/certificates experiments/instances
cp "$tmp"/certificates/* experiments/certificates/
cp "$tmp"/instances/* experiments/instances/
rm -rf "$tmp"
echo "artifacts for $TAG fetched and checked against experiments/ARTIFACTS.sha256"
