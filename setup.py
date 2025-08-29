#!/usr/bin/env python3
"""
Setup script for BELT (BalazsLab ERRSeq Long-Read Toolkit)
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "BELT - BalazsLab ERRSeq Long-Read Toolkit"

setup(
    name="belt",
    version="0.1.0",
    description="BalazsLab ERRSeq Long-Read Toolkit",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    author="BalazsLab",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/belt",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click>=8.0.0",
        "pyyaml>=6.0",
    ],
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'belt=__main__:main',
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics, sequencing, UMI, long-read",
)
