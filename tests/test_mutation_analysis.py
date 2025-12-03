#!/usr/bin/env python3
"""
Test script for BELT mutation analysis functionality
"""

import os
import sys
from pathlib import Path

# Add scripts directory to path
scripts_dir = Path(__file__).parent / "scripts"
sys.path.insert(0, str(scripts_dir))

def test_mutation_analysis():
    """Test the mutation analysis functionality"""
    print("Testing BELT Mutation Analysis...")
    
    # Check if we can import the module
    try:
        from mutation_analysis import analyze_mutations
        print("✅ Mutation analysis module imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import mutation analysis module: {e}")
        return False
    
    # Check if we can import the web app
    try:
        import mutation_visualizer
        print("✅ Mutation visualizer module imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import mutation visualizer module: {e}")
        return False
    
    print("✅ All modules imported successfully")
    return True

if __name__ == "__main__":
    success = test_mutation_analysis()
    if success:
        print("\n🎉 Mutation analysis functionality is ready!")
        print("\nTo use the mutation analysis:")
        print("1. Run: belt mutation-analyze --consensus-dir <dir> --reference-dir <dir> --output-dir <dir>")
        print("2. Run: python mutation_visualizer.py")
        print("3. Open browser to http://localhost:8050")
    else:
        print("\n❌ Mutation analysis functionality has issues")
        sys.exit(1)
