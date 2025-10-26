#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

DATA="./data"
BIN_DIR="./bin"
CODON_BIN="$BIN_DIR/align_codon"

# --- helpers ---------------------------------------------------------------

time_ms() {
  python - "$@" <<'PY'
import subprocess, sys, time
t=time.time()
subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
print(int((time.time()-t)*1000))
PY
}

print_row() { printf "%-18s%-12s%sms\n" "$1" "$2" "$3"; }

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "Missing required file: $1" 1>&2
    exit 1
  fi
}

ensure_libomp_next_to_binary() {
  # If already present, we're good.
  if [[ -f "$BIN_DIR/libomp.dylib" ]]; then
    return
  fi

  # Try Homebrew first (Apple Silicon: /opt/homebrew, Intel: /usr/local)
  if command -v brew >/dev/null 2>&1; then
    OMP="$(brew --prefix libomp 2>/dev/null || true)/lib/libomp.dylib"
    if [[ -f "$OMP" ]]; then
      ln -sf "$OMP" "$BIN_DIR/libomp.dylib"
      return
    fi
  fi

  # Common fallback locations
  for P in \
    /opt/homebrew/opt/libomp/lib/libomp.dylib \
    /usr/local/opt/libomp/lib/libomp.dylib
  do
    if [[ -f "$P" ]]; then
      ln -sf "$P" "$BIN_DIR/libomp.dylib"
      return
    fi
  done

  echo "Warning: libomp.dylib not found. If the binary fails to run, install libomp (e.g. 'brew install libomp')." 1>&2
}

# --- sanity checks ---------------------------------------------------------

# required data (you said these are already in repo)
require_file "$DATA/MT-human.fa"
require_file "$DATA/MT-orang.fa"
require_file "$DATA/q1.fa"
require_file "$DATA/t1.fa"
require_file "$DATA/q2.fa"
require_file "$DATA/t2.fa"

# --- build codon binary if missing (hard-require compiled binary) ----------

if [[ ! -x "$CODON_BIN" ]]; then
  if ! command -v codon >/dev/null 2>&1; then
    echo "Error: Codon not found, and required binary '$CODON_BIN' is missing." 1>&2
    echo "Install Codon and build: codon build -release -o $CODON_BIN codon_impl/align_codon.py" 1>&2
    exit 1
  fi
  mkdir -p "$BIN_DIR"
  # Try to build once
  codon build -release -o "$CODON_BIN" codon_impl/align_codon.py
fi

# Ensure libomp.dylib sits next to the binary for macOS loaders
ensure_libomp_next_to_binary
# (extra safety) make sure the loader also sees bin/ first
export DYLD_LIBRARY_PATH="$PWD/$BIN_DIR:${DYLD_LIBRARY_PATH:-}"

# --- run table -------------------------------------------------------------

echo "Method            Language    Runtime"
echo "--------------------------------------"

# 1) global-mt_human (python)
ms=$(time_ms python -m py_impl.align_cli \
  --method global \
  --query "$DATA/MT-human.fa" \
  --target "$DATA/MT-orang.fa")
print_row "global-mt_human" "python" "$ms"

# 2) global-q1 (python)
ms=$(time_ms python -m py_impl.align_cli \
  --method global \
  --query "$DATA/q1.fa" \
  --target "$DATA/t1.fa")
print_row "global-q1" "python" "$ms"

# 3) local-q2 (Codon ONLY; compiled binary)
ms=$(time_ms "$CODON_BIN" \
  --method local \
  --query "$DATA/q2.fa" \
  --target "$DATA/t2.fa")
print_row "local-q2" "codon" "$ms"
