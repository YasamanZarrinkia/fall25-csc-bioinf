#!/usr/bin/env bash
set -euo pipefail

echo "🚀 Starting Week 2 Evaluation – Codon & Python (running inside week2/)"

# Always execute from week2 so relative paths (and artifacts) live here
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CODE_DIR="code"

echo "📁 Checking required files..."
required_files=(
  "$CODE_DIR/decomposer.py"
  "$CODE_DIR/main.py"
  "$CODE_DIR/test_codon.py"
  "$CODE_DIR/test_python.py"
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
echo "🧪 Running Codon tests ($CODE_DIR/test_codon.py)…"
PYTHONPATH="${PYTHONPATH:-}:$CODE_DIR" codon run "$CODE_DIR/test_codon.py"
echo "✅ Codon tests passed"

echo ""
echo "🐍 Running Python tests ($CODE_DIR/test_python.py)…"
PYTHONPATH="${PYTHONPATH:-}:$CODE_DIR" python "$CODE_DIR/test_python.py"
echo "✅ Python tests passed"

echo ""
echo "🧹 Cleaning temporary artifacts (motif map files, if any)…"
rm -f ./*_motif_map.txt || true

echo ""
echo "🎉 WEEK 2 EVALUATION COMPLETE"
echo "📊 SUMMARY:"
echo "   ✅ All required files present"
echo "   ✅ Codon test script passed"
echo "   ✅ Python test script passed"
