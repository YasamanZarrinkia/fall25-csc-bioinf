#!/bin/bash

set -e  # Exit on any error

echo "🚀 Starting Week 2 Evaluation - Codon TRViz Implementation"

# Change to code directory
cd code

echo "📁 Checking required files..."
required_files=("decomposer.py" "main.py" "test_codon.py" "test_python.py")
missing_files=0
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ Found: $file"
    else
        echo "❌ Missing: $file"
        missing_files=1
    fi
done

if [ $missing_files -ne 0 ]; then
    echo "❌ Missing required files. Evaluation failed."
    exit 1
fi

echo ""
echo "🧪 Running Codon tests..."
if codon run test_codon.py; then
    echo "✅ Codon tests passed"
else
    echo "❌ Codon tests failed"
    exit 1
fi

echo ""
echo "🐍 Running Python tests..."
if python test_python.py; then
    echo "✅ Python tests passed"
else
    echo "❌ Python tests failed"
    exit 1
fi

echo ""
echo "🔧 Testing basic decomposition functionality..."
python -c "
from decomposer import Decomposer
print('Testing perfect repeats...')
d = Decomposer('DP')
result = d.decompose('ACTGACTGACTG', ['ACTG'])
print(f'Result: {result}')
assert result == ['ACTG', 'ACTG', 'ACTG'], f'Expected different decomposition: {result}'

print('Testing imperfect repeats...')
result2 = d.decompose('ACTGACTTACTG', ['ACTG'])
print(f'Result: {result2}')
assert result2 == ['ACTG', 'ACTT', 'ACTG'], f'Expected different decomposition: {result2}'
print('✅ Basic decomposition tests passed')
"

echo ""
echo "🔄 Testing integration pipeline..."
python -c "
from main import TandemRepeatVizWorker
worker = TandemRepeatVizWorker()
print('Testing small integration...')
result = worker.decomposer.decompose('ACTACT', ['ACT'])
print(f'Integration test result: {result}')
assert result == ['ACT', 'ACT'], f'Integration test failed: {result}'

print('Testing refinement...')
decomposed_trs = [['ACT', 'ACT'], ['ACT', 'ACC', 'ACT']]
refined = worker.decomposer.refine(decomposed_trs)
print(f'Refinement result: {refined}')
assert len(refined) == 2, 'Refinement failed'
print('✅ Integration tests passed')
"

echo ""
echo "🎉 WEEK 2 EVALUATION COMPLETE!"
echo "📊 SUMMARY:"
echo "   ✅ All required files present"
echo "   ✅ Codon test suite passed"
echo "   ✅ Python test suite passed" 
echo "   ✅ Basic decomposition functionality verified"
echo "   ✅ Integration pipeline working"
echo ""
echo "🚀 Codon TRViz implementation is fully validated and ready for submission!"
