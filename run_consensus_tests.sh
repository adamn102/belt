#!/bin/bash
# Script to run consensus tests specifically

echo "================================================"
echo "Running BELT Consensus Tests"
echo "================================================"

# Navigate to the belt directory
cd /Users/anitido/Documents/GitHub/belt/belt

echo ""
echo "Activating belt-pipeline-medaka environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate belt-pipeline-medaka

echo ""
echo "Running consensus tests only..."
echo ""

# Run only the consensus tests
pytest tests/test_consensus.py -v --tb=short -s

TEST_EXIT_CODE=$?

echo ""
echo "================================================"
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "✅ All consensus tests passed!"
else
    echo "❌ Some consensus tests failed (exit code: $TEST_EXIT_CODE)"
fi
echo "================================================"

exit $TEST_EXIT_CODE





