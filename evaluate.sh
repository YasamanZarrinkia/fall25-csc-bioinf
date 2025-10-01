#!/usr/bin/env bash
set -euo pipefail

echo "🚀 Starting Week 2 Evaluation – Codon & Python"

# Ensure repo layout
if [ ! -d "code" ]; then
  echo "❌ ./code directory not found (run from repo root)."
  exit 1
fi

echo "📁 Checking required files..."
required_files=(
  "code/decomposer.py"
  "code/main.py"
  "code/test_codon.py"
  "code/test_python.py"
)
missing=0
for f in "${required_files[@]}"; do
  if [ -f "$f" ]; then
    echo "✅ Found: $f"
  else
    echo "❌ Missing: $f"
    missing=1
  fi
done
if [ $missing -ne 0 ]; then
  echo "❌ Missing required files. Evaluation failed."
  exit 1
fi

echo ""
echo "ℹ️ Tool versions"
python --version || true
pip --version || true
codon --version || true

echo ""
echo "🧪 Running Codon tests (code/test_codon.py)…"
# Ensure local imports resolve (your test_codon imports from same dir)
PYTHONPATH="${PYTHONPATH:-}:code" codon run code/test_codon.py
echo "✅ Codon tests passed"

echo ""
echo "🐍 Running Python tests (code/test_python.py)…"
# Your Python test already appends its dirname to sys.path, but we also set PYTHONPATH for safety
PYTHONPATH="${PYTHONPATH:-}:code" python code/test_python.py
echo "✅ Python tests passed"

echo ""
echo "🧹 Cleaning temporary artifacts (motif map files, if any)…"
rm -f ./test_motif_map.txt ./*_motif_map.txt || true

echo ""
echo "🎉 WEEK 2 EVALUATION COMPLETE"
echo "📊 SUMMARY:"
echo "   ✅ All required files present"
echo "   ✅ Codon test script passed"
echo "   ✅ Python test script passed"
