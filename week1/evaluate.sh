#!/usr/bin/env bash
set -euo pipefail

### ── Locate repo root (works from anywhere) ────────────────────────────────────
if ROOT="$(git rev-parse --show-toplevel 2>/dev/null)"; then
  :
else
  SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
  ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

### ── Config (override via env) ────────────────────────────────────────────────
PY_DIR="${PY_DIR:-$ROOT/week1/code/genome-assembly}"
CODON_DIR="${CODON_DIR:-$ROOT/week1/code/genome-assembly-codon}"
DATA_BASE="${DATA_BASE:-$ROOT/week1/data}"
OUT_DIR="${OUT_DIR:-$ROOT/tmp}"
DATASETS=(${DATASETS:-data1 data2 data3 data4})

# Timeouts (seconds)
PY_TIMEOUT_SEC="${PY_TIMEOUT_SEC:-1800}"        # 30m
CODON_TIMEOUT_SEC="${CODON_TIMEOUT_SEC:-1200}"  # 20m

# Docker controls
USE_DOCKER_PY="${USE_DOCKER_PY:-auto}"               # auto|true|false
USE_DOCKER_CODON="${USE_DOCKER_CODON:-auto}"         # auto|true|false
DOCKER_PLATFORM_PY="${DOCKER_PLATFORM_PY:-}"         # e.g. linux/amd64
DOCKER_PLATFORM_CODON="${DOCKER_PLATFORM_CODON:-linux/amd64}"  # safer on Apple Silicon
DOCKER_CPUS="${DOCKER_CPUS:-4}"
DOCKER_MEM="${DOCKER_MEM:-8g}"

mkdir -p "$OUT_DIR"

### ── Helpers ─────────────────────────────────────────────────────────────────
have() { command -v "$1" >/dev/null 2>&1; }
fmt_time(){ s=${1%.*}; h=$((s/3600)); m=$(((s%3600)/60)); s=$((s%60)); printf "%d:%02d:%02d" "$h" "$m" "$s"; }
run_with_time(){ # $1=cmd $2=stdout $3=stderr
  local t; t="$(mktemp)"
  if /usr/bin/time -p -o "$t" bash -lc "$1" >"${2:-/dev/null}" 2>"${3:-/dev/null}"; then :; else :; fi
  awk '/^real /{print $2}' "$t" 2>/dev/null || true
  rm -f "$t"
}

n50_of_fasta(){ python3 - "$1" <<'PY' 2>/dev/null || true
import sys,re,os
p=sys.argv[1]
if not os.path.isfile(p): print("NA"); raise SystemExit
lens=[]; cur=None
with open(p,'r',errors='ignore') as f:
  for line in f:
    s=line.strip()
    if not s: continue
    if s.startswith('>'):
      if cur: lens.append(cur)
      cur=0
    else:
      if re.fullmatch(r'[ACGTNacgtn]+', s): cur=(cur or 0)+len(s)
if cur: lens.append(cur)
if not lens: print("NA"); raise SystemExit
lens.sort(reverse=True); tot=sum(lens); half=tot/2; acc=0
for L in lens:
  acc+=L
  if acc>=half: print(L); break
PY
}

prepare_python_local(){
  if have python3; then
    python3 -m pip install --upgrade pip >/dev/null 2>&1 || true
    # matplotlib only if your Python code imports it
    python3 -m pip install --no-cache-dir matplotlib >/dev/null 2>&1 || true
  fi
}

prepare_codon_local(){
  if ! have codon && have curl; then
    /bin/bash -c "$(curl -fsSL https://exaloop.io/install.sh)" || true
    export PATH="$HOME/.codon/bin:$PATH"
  fi
}

### ── Python runners ───────────────────────────────────────────────────────────
run_python_local(){
  local ds="$1" d="$DATA_BASE/$ds" of="$OUT_DIR/${ds}.py.stdout" ef="$OUT_DIR/${ds}.py.stderr"
  [ -d "$d" ] || { echo "-,NA"; return; }
  rm -f "$d/contig.fasta" "$d/contig_py.fasta"
  (ulimit -s 8192000 || true) >/dev/null 2>&1
  local secs; secs="$(run_with_time "timeout $PY_TIMEOUT_SEC python3 \"$PY_DIR/main.py\" \"$d\"" "$of" "$ef" || true)"
  [ -f "$d/contig.fasta" ] && mv -f "$d/contig.fasta" "$d/contig_py.fasta" || true
  local pretty="-"; [ -n "${secs:-}" ] && pretty="$(fmt_time "$secs")"
  echo "${pretty},$(n50_of_fasta "$d/contig_py.fasta")"
}

run_python_docker(){
  local ds="$1" d="$DATA_BASE/$ds" of="$OUT_DIR/${ds}.py.stdout" ef="$OUT_DIR/${ds}.py.stderr"
  [ -d "$d" ] || { echo "-,NA"; return; }
  rm -f "$d/contig.fasta" "$d/contig_py.fasta"
  local plat=""; [ -n "$DOCKER_PLATFORM_PY" ] && plat="--platform=$DOCKER_PLATFORM_PY"
  local cmd="pip install -q matplotlib && MPLBACKEND=Agg timeout $PY_TIMEOUT_SEC python3 main.py \"$d\""
  local secs; secs="$(run_with_time "docker run --rm $plat -v \"$ROOT\":\"$ROOT\" -w \"$PY_DIR\" python:3.11-slim bash -lc '$cmd'" "$of" "$ef" || true)"
  [ -f "$d/contig.fasta" ] && mv -f "$d/contig.fasta" "$d/contig_py.fasta" || true
  local pretty="-"; [ -n "${secs:-}" ] && pretty="$(fmt_time "$secs")"
  echo "${pretty},$(n50_of_fasta "$d/contig_py.fasta")"
}

run_python_dataset(){
  local ds="$1"
  case "$USE_DOCKER_PY" in
    true)  run_python_docker "$ds" ;;
    false) run_python_local  "$ds" ;;
    *)     if have python3; then run_python_local "$ds"
           elif have docker;  then run_python_docker "$ds"
           else echo "-,NA"; fi ;;
  esac
}

### ── Codon runners ────────────────────────────────────────────────────────────
run_codon_local(){
  local ds="$1" d="$DATA_BASE/$ds" of="$OUT_DIR/${ds}.codon.stdout" ef="$OUT_DIR/${ds}.codon.stderr"
  [ -d "$d" ] || { echo "-,NA"; return; }
  rm -f "$d/contig.fasta" "$d/contig_codon.fasta"
  local secs; secs="$(run_with_time "timeout $CODON_TIMEOUT_SEC codon run -release \"$CODON_DIR/main.py\" \"$d\"" "$of" "$ef" || true)"
  if [ -z "${secs:-}" ]; then
    secs="$(run_with_time "timeout $CODON_TIMEOUT_SEC codon run \"$CODON_DIR/main.py\" \"$d\"" "$of" "$ef" || true)"
  fi
  [ -f "$d/contig.fasta" ] && mv -f "$d/contig.fasta" "$d/contig_codon.fasta" || true
  local pretty="-"; [ -n "${secs:-}" ] && pretty="$(fmt_time "$secs")"
  echo "${pretty},$(n50_of_fasta "$d/contig_codon.fasta")"
}

run_codon_docker(){
  local ds="$1" d="$DATA_BASE/$ds" of="$OUT_DIR/${ds}.codon.stdout" ef="$OUT_DIR/${ds}.codon.stderr"
  [ -d "$d" ] || { echo "-,NA"; return; }
  rm -f "$d/contig.fasta" "$d/contig_codon.fasta"
  local plat="--platform=${DOCKER_PLATFORM_CODON}"
  local res="--cpus=${DOCKER_CPUS} --memory=${DOCKER_MEM}"
  local cmd="timeout $CODON_TIMEOUT_SEC codon run -release main.py \"$d\""
  local secs; secs="$(run_with_time "docker run --rm $plat $res -v \"$ROOT\":\"$ROOT\" -w \"$CODON_DIR\" exaloop/codon:latest bash -lc '$cmd'" "$of" "$ef" || true)"
  if [ -z "${secs:-}" ]; then
    cmd="timeout $CODON_TIMEOUT_SEC codon run main.py \"$d\""
    secs="$(run_with_time "docker run --rm $plat $res -v \"$ROOT\":\"$ROOT\" -w \"$CODON_DIR\" exaloop/codon:latest bash -lc '$cmd'" "$of" "$ef" || true)"
  fi
  [ -f "$d/contig.fasta" ] && mv -f "$d/contig.fasta" "$d/contig_codon.fasta" || true
  local pretty="-"; [ -n "${secs:-}" ] && pretty="$(fmt_time "$secs")"
  echo "${pretty},$(n50_of_fasta "$d/contig_codon.fasta")"
}

run_codon_dataset(){
  local ds="$1"
  case "$USE_DOCKER_CODON" in
    true)  run_codon_docker "$ds" ;;
    false) run_codon_local  "$ds" ;;
    *)     if have codon; then run_codon_local "$ds"
           elif have docker; then run_codon_docker "$ds"
           else echo "-,NA"; fi ;;
  esac
}

### ── Bootstrap locals ─────────────────────────────────────────────────────────
prepare_python_local
prepare_codon_local

### ── Header + runs ────────────────────────────────────────────────────────────
printf '%s\n' "Dataset	Language	Runtime	N50"
printf '%s\n' "-------------------------------------------------------------------------------------------------------"

for ds in "${DATASETS[@]}"; do
  IFS=, read -r py_rt py_n50 <<<"$(run_python_dataset "$ds")"
  printf '%s\t%s\t\t%s\t%s\n' "$ds" "python" "${py_rt:-"-"}" "${py_n50:-"NA"}"

  IFS=, read -r co_rt co_n50 <<<"$(run_codon_dataset "$ds")"
  printf '%s\t%s\t\t%s\t%s\n' "$ds" "codon"  "${co_rt:-"-"}" "${co_n50:-"NA"}"
done

exit 0
