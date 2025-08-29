#!/bin/bash

# BELT Development Environment Setup Script

echo "Setting up BELT development environment..."

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv .venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source .venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing requirements..."
pip install -r requirements.txt

# Install the package in development mode
echo "Installing BELT in development mode..."
pip install -e .

echo ""
echo "✅ BELT development environment setup complete!"
echo ""
echo "To use BELT:"
echo "1. Activate the virtual environment: source .venv/bin/activate"
echo "2. Run belt commands: belt --help"
echo ""
echo "Available commands:"
echo "  belt --help                    # Show all commands"
echo "  belt extract --help            # Show extract command help"
echo "  belt umi-map --help            # Show umi-map command help"
echo "  belt umi-group --help          # Show umi-group command help"
echo "  belt whitelist --help          # Show whitelist command help"
echo "  belt run-pipeline --help       # Show run-pipeline command help"
