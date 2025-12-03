#!/bin/bash
# Script to run belt tests in the correct conda environment

echo "================================================"
echo "Running BELT Tests"
echo "================================================"

# Check if conda environment exists
if ! conda env list | grep -q "belt-pipeline-medaka"; then
    echo "Error: belt-pipeline-medaka conda environment not found"
    echo "Please create it first or activate the correct environment"
    exit 1
fi

# Navigate to the belt directory
cd /Users/anitido/Documents/GitHub/belt/belt

echo ""
echo "Activating belt-pipeline-medaka environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate belt-pipeline-medaka

echo ""
echo "Running pytest..."
echo ""

# Run all tests
pytest tests/ -v --tb=short

TEST_EXIT_CODE=$?

echo ""
echo "================================================"
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "✅ All tests passed!"
else
    echo "❌ Some tests failed (exit code: $TEST_EXIT_CODE)"
fi
echo "================================================"

exit $TEST_EXIT_CODE





