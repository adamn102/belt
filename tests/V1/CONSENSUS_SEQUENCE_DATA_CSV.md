# Consensus Sequence Data CSV Generation

## Overview

Added functionality to `mutation_summary.py` to generate a detailed CSV file containing consensus sequence data for each sample, including UMI IDs, sequence lengths, and mutation counts.

## New Functionality

### Function Added

**`generate_consensus_sequence_data(sample_results: List[Dict], output_dir: Path)`**

Generates a CSV file with detailed information about each consensus sequence.

### CSV Output

**File**: `consensus_sequence_data.csv`

**Location**: `data/10-mutation_analysis/consensus_sequence_data.csv`

**Columns**:
- `Sample`: Sample name
- `UMI_ID`: Consensus sequence identifier (UMI ID)
- `Total_Bases`: Length of the consensus sequence
- `Total_Mutations`: Total number of mutations in this consensus
- `Mutation_Rate`: Mutation rate (mutations per base)
- `Substitutions`: Number of substitution mutations
- `Insertions`: Number of insertion mutations
- `Deletions`: Number of deletion mutations

### Example CSV Output

```csv
Sample,UMI_ID,Total_Bases,Total_Mutations,Mutation_Rate,Substitutions,Insertions,Deletions
Sample_A,UMI_001,1000,5,0.005,2,1,2
Sample_A,UMI_002,950,2,0.0021,1,1,0
Sample_B,UMI_003,1200,8,0.0067,3,2,3
```

## Integration

### Pipeline Integration

The new CSV file is automatically generated as part of Step 10 (mutation analysis) in the Belt Pipeline.

**Updated Snakefile**:
- Added `consensus_sequence_data.csv` to `list_mutation_analysis_outputs()`
- File is now included in pipeline dependencies and outputs

### Function Integration

The `generate_consensus_sequence_data()` function is called automatically within the `mutation_summary()` function, so no changes are needed to existing pipeline commands.

## Usage

### Automatic Generation

The CSV is generated automatically when running the mutation analysis step:

```bash
# Via pipeline
snakemake data/10-mutation_analysis/consensus_sequence_data.csv

# Via belt command
belt mutation-summary --output-dir data/10-mutation_analysis
```

### Manual Generation

You can also generate the CSV manually:

```python
from mutation_summary import generate_consensus_sequence_data
from pathlib import Path

# Load sample results (from JSON files)
sample_results = [...]  # List of sample result dictionaries

# Generate CSV
generate_consensus_sequence_data(sample_results, Path("output_dir"))
```

## Data Source

The function extracts data from the `per_read_analysis` field in each sample result, which contains:

```python
{
    'consensus_id': 'UMI_001',
    'sequence_length': 1000,
    'mutation_count': 5,
    'mutation_rate': 0.005,
    'mutations': [
        {'type': 'substitution'},
        {'type': 'substitution'},
        {'type': 'insertion'},
        {'type': 'deletion'},
        {'type': 'deletion'}
    ]
}
```

## Error Handling

The function handles various edge cases:

1. **Empty sample results**: Creates empty CSV with proper headers
2. **Missing per_read_analysis**: Skips samples without consensus data
3. **Missing mutations data**: Defaults mutation type counts to 0
4. **No consensus sequences**: Creates empty CSV file

## Test Files

### Test Suite

**`tests/test_consensus_sequence_data.py`**
- Comprehensive unittest suite
- Tests normal operation, edge cases, and error handling
- Verifies CSV structure and content

### Verification Script

**`tests/verify_consensus_csv.py`**
- Simple verification script
- Creates mock data and tests CSV generation
- Provides clear pass/fail output

### Running Tests

```bash
cd /Users/anitido/Documents/GitHub/belt/belt

# Run comprehensive test suite
python3 -m unittest tests.test_consensus_sequence_data -v

# Run quick verification
python3 tests/verify_consensus_csv.py
```

## Benefits

1. **Detailed Analysis**: Per-consensus sequence data for detailed analysis
2. **UMI Tracking**: Links mutations to specific UMI groups
3. **Mutation Breakdown**: Separate counts for substitutions, insertions, deletions
4. **Easy Integration**: Works with existing analysis tools (pandas, R, etc.)
5. **Pipeline Integration**: Automatically generated as part of standard workflow

## Use Cases

1. **Quality Control**: Identify consensus sequences with high mutation rates
2. **UMI Analysis**: Track mutation patterns across UMI groups
3. **Sample Comparison**: Compare mutation rates between samples
4. **Downstream Analysis**: Import into statistical analysis tools
5. **Visualization**: Create plots showing mutation distribution

## Files Modified

### Core Code
- `scripts/mutation_summary.py`: Added `generate_consensus_sequence_data()` function

### Pipeline
- `Belt Pipeline/Snakefile`: Updated `list_mutation_analysis_outputs()` to include new CSV

### Test Files
- `tests/test_consensus_sequence_data.py`: Comprehensive test suite
- `tests/verify_consensus_csv.py`: Simple verification script

### Documentation
- `tests/CONSENSUS_SEQUENCE_DATA_CSV.md`: This documentation file

## Date

Implemented: October 28, 2025
