#!/usr/bin/env python3
"""
BELT CLI module
"""

import sys
import os

# Add the scripts directory to the path
scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts')
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

# Import the main function from belt.py
from belt import main

if __name__ == "__main__":
    main()
