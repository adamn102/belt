# BELT Development Guide

This guide explains how to set up and use the BELT (BalazsLab ERRSeq Long-Read Toolkit) development environment.

## Quick Start

### 1. Clone the Repository
```bash
git clone <your-repo-url>
cd belt
```

### 2. Set Up Development Environment
```bash
# Run the automated setup script
./install.sh
```

### 3. Activate the Virtual Environment
```bash
source .venv/bin/activate
```

### 4. Use BELT
```bash
# Run belt commands
./belt.py --help
# or
./belt-cmd --help
```

## Manual Setup (Alternative)

If you prefer to set up manually:

### 1. Create Virtual Environment
```bash
python3 -m venv .venv
source .venv/bin/activate  # On Unix/macOS
# or
.venv\Scripts\activate     # On Windows
```

### 2. Install Dependencies
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Install BELT in Development Mode
```bash
pip install -e .
```

## Using BELT Commands

Once set up, you can use the following commands:

### Main Help
```bash
./belt.py --help
```

### Available Commands

#### 1. Extract Command
```bash
./belt.py extract --help
./belt.py extract --reference <ref.fasta> --alignment <input.bam> --output <output.bam>
```

#### 2. UMI Map Command
```bash
./belt.py umi-map --help
./belt.py umi-map --whitelist <whitelist.csv> --reference <ref.fasta> --alignment <input.bam> --motif <motif> --output_dir <output_dir>
```

#### 3. UMI Group Command
```bash
./belt.py umi-group --help
./belt.py umi-group --reference <ref.fasta> --alignment <input.bam> --corrected_umis <umis.csv> --output_dir <output_dir> --grouping 5pUMI
```

#### 4. Whitelist Command
```bash
./belt.py whitelist --help
./belt.py whitelist --reference <ref.fasta> --alignment <input.bam> --motif <motif> --output_dir <output_dir>
```

#### 5. Consensus Command
```bash
./belt.py consensus --help
./belt.py consensus --grouped-bam <grouped.bam> --reference <ref.fasta> --output-dir <output_dir> --sample <sample_name>
```

#### 6. Run Pipeline Command
```bash
./belt.py run-pipeline --help
./belt.py run-pipeline --config <config.yaml>
```

## Development Workflow

### 1. Activate Environment
Always activate the virtual environment before working:
```bash
source .venv/bin/activate
```

### 2. Make Changes
Edit the files in the `scripts/` directory as needed.

### 3. Test Changes
```bash
# Test a specific command
./belt.py extract --help

# Test with sample data
./belt.py extract --reference test-data/test_ref.fasta --alignment test-data/test_sorted.bam --output test_output.bam
```

### 4. Install Updates
After making changes to the code:
```bash
pip install -e .
```

## Project Structure

```
belt/
├── belt.py              # Main entry point
├── belt-cmd             # Command wrapper script
├── install.sh           # Setup script
├── requirements.txt     # Python dependencies
├── setup.py            # Package configuration
├── scripts/            # Main BELT scripts
│   ├── belt.py         # CLI implementation
│   ├── extract.py      # Extract functionality
│   ├── umi_map.py      # UMI mapping
│   ├── umi_group.py    # UMI grouping
│   ├── whitelist.py    # Whitelist generation
│   ├── consensus.py    # Consensus calling
│   └── utils.py        # Utility functions
├── test-data/          # Test data files
├── tests/              # Test files
└── .venv/              # Virtual environment (created by setup)
```

## Troubleshooting

### Command Not Found
If you get "command not found" for belt:
1. Make sure the virtual environment is activated
2. Try running `./belt.py` instead of `belt`
3. Check that `install.sh` completed successfully

### Import Errors
If you get import errors:
1. Make sure all dependencies are installed: `pip install -r requirements.txt`
2. Reinstall the package: `pip install -e .`
3. Check that the virtual environment is activated

### Permission Errors
If you get permission errors:
```bash
chmod +x belt.py
chmod +x belt-cmd
chmod +x install.sh
```

## Adding New Commands

To add a new command to BELT:

1. Add the command function to `scripts/belt.py`:
```python
@main.command()
@click.option('--input', required=True, help="Input file")
def new_command(input):
    """Description of new command"""
    # Your command logic here
    pass
```

2. Test the new command:
```bash
./belt.py new-command --help
```

## Testing

Run tests to ensure everything works:
```bash
# Run all tests
python -m pytest tests/

# Run specific test file
python -m pytest tests/test_belt_cli_simple.py

# Run with coverage
python -m pytest tests/ --cov=scripts
```

## Contributing

1. Create a feature branch
2. Make your changes
3. Test thoroughly
4. Submit a pull request

## Environment Variables

You can set these environment variables for development:

- `BELT_DEBUG=1` - Enable debug output
- `BELT_LOG_LEVEL=DEBUG` - Set log level

## IDE Setup

### VS Code
1. Open the project folder
2. Select the Python interpreter from `.venv/bin/python`
3. Install the Python extension

### PyCharm
1. Open the project
2. Go to Settings > Project > Python Interpreter
3. Add the interpreter from `.venv/bin/python`

## Deployment

To deploy BELT to other systems:

1. Package the application:
```bash
python setup.py sdist bdist_wheel
```

2. Install on target system:
```bash
pip install dist/belt-0.1.0.tar.gz
```

## Support

For issues and questions:
1. Check this development guide
2. Look at the test files for examples
3. Create an issue in the repository
