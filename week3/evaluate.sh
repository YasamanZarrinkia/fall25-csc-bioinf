#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pushd "$SCRIPT_DIR" >/dev/null

PY_OUT="$(python3 test_py.py)"
CD_OUT="$(codon run -release test_codon.py)"

echo "Language    Runtime"
echo "-------------------"
echo "$PY_OUT" | awk '/^python[[:space:]]/{print;exit}'
echo "$CD_OUT" | awk '/^codon[[:space:]]/{print;exit}'

popd >/dev/null
