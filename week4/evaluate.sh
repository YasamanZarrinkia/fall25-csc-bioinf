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

# macOS-only helper (no-op on Linux)
ensure_libomp_next_to_binary_macos() {
  [[ "$(uname -s)" != "Darwin" ]] && return 0
  [[ -f "$BIN_DIR/libomp.dylib" ]] && return 0

  if command -v brew >/dev/null 2>&1; then
    OMP="$(brew --prefix libomp 2>/dev/null || true)/lib/libomp.dylib"
    if [[ -f "$OMP" ]]; then
      ln -sf "$OMP" "$BIN_DIR/libomp.dylib"
      return 0
    fi
  fi
  for P in /opt/homebrew/opt/libomp/lib/libomp.dylib /usr/local/opt/libomp/lib/libomp.dylib; do
    [[ -f "$P" ]] && { ln -sf "$P" "$BIN_DIR/libomp.dylib"; return 0; }
  done
  echo "Warning: libomp.dylib not found on macOS. Try 'brew install libomp'." 1>&2
}

# --- sanity checks ---------------------------------------------------------

require_file "$DATA/MT-human.fa"
require_file "$DATA/MT-orang.fa"
require_file "$DATA/q1.fa"
require_file "$DATA/t1.fa"
require_file "$DATA/q2.fa"
require_file "$DATA/t2.fa"

# --- build codon binary if missing (build per-OS to avoid Exec format errors)

if [[ "$(uname -s)" == "Linux" && -x "$CODON_BIN" ]]; then
  rm -f "$CODON_BIN"
fi

if [[ ! -x "$CODON_BIN" ]]; then
  if ! command -v codon >/dev/null 2>&1; then
    echo "Error: Codon not found, and required binary '$CODON_BIN' is missing." 1>&2
    echo "Install Codon and build: codon build -release -o $CODON_BIN codon_impl/align_codon.py" 1>&2
    exit 1
  fi
  mkdir -p "$BIN_DIR"
  codon build -release -o "$CODON_BIN" codon_impl/align_codon.py
fi

# macOS needs libomp alongside the binary; Linux typically finds libomp.so via system paths
ensure_libomp_next_to_binary_macos

# Export platform-specific loader path only when needed
if [[ "$(uname -s)" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="$PWD/$BIN_DIR:${DYLD_LIBRARY_PATH:-}"
else
  # Usually not necessary on Ubuntu runners, but harmless
  export LD_LIBRARY_PATH="$PWD/$BIN_DIR:${LD_LIBRARY_PATH:-}"
fi

# --- run table -------------------------------------------------------------

echo "Method            Language    Runtime"
echo "--------------------------------------"

ms=$(time_ms python -m py_impl.align_cli \
  --method global \
  --query "$DATA/MT-human.fa" \
  --target "$DATA/MT-orang.fa")
print_row "global-mt_human" "python" "$ms"

ms=$(time_ms python -m py_impl.align_cli \
  --method global \
  --query "$DATA/q1.fa" \
  --target "$DATA/t1.fa")
print_row "global-q1" "python" "$ms"

ms=$(time_ms "$CODON_BIN" \
  --method local \
  --query "$DATA/q2.fa" \
  --target "$DATA/t2.fa")
print_row "local-q2" "codon" "$ms"
