#!/bin/bash
# Helper script to run UMI grouping tests and save output for visual inspection

echo "Running UMI Grouping Tests with Visual Sequence Output"
echo "========================================================"
echo ""

cd "$(dirname "$0")/.."

# Run all tests with verbose output
echo "Running all UMI grouping tests..."
pytest tests/test_umi_group.py -v -s --tb=short | tee tests/umi_group_test_output.txt

echo ""
echo "========================================================"
echo "Test output saved to: tests/umi_group_test_output.txt"
echo ""
echo "To run individual tests:"
echo "  pytest tests/test_umi_group.py::test_umi_grouping_with_fastq_inspection -v -s"
echo "  pytest tests/test_umi_group.py::test_umi_grouping_different_columns -v -s"
echo "  pytest tests/test_umi_group.py::test_umi_grouping_fast -v -s"
echo "  pytest tests/test_umi_group.py::test_umi_grouping_fast_with_fastq -v -s"

