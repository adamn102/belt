#!/usr/bin/env python3
"""
BELT UMI Meta-Analysis Module
Analyzes UMI patterns and statistics from UMI mapping data
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional
import logging
from collections import Counter
import re

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def analyze_umi_patterns(umi_file: Path) -> Dict:
    """
    Analyze UMI patterns for a single sample
    
    Args:
        umi_file: Path to UMI mapping CSV file
    
    Returns:
        Dictionary containing UMI analysis results
    """
    logger.info(f"Analyzing UMI patterns for: {umi_file}")
    
    try:
        # Read the UMI mapping file
        df = pd.read_csv(umi_file)
        
        # Extract UMI sequences
        umi_sequences = []
        if 'UMI1' in df.columns:
            umi_sequences.extend(df['UMI1'].dropna().tolist())
        if 'UMI2' in df.columns:
            umi_sequences.extend(df['UMI2'].dropna().tolist())
        
        if not umi_sequences:
            logger.warning(f"No UMI sequences found in {umi_file}")
            return {
                'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
                'total_umis': 0,
                'unique_umis': 0,
                'umi_length_stats': {},
                'umi_composition': {},
                'umi_diversity': 0,
                'read_counts': {}
            }
        
        # Basic statistics
        total_umis = len(umi_sequences)
        unique_umis = len(set(umi_sequences))
        
        # UMI length statistics
        umi_lengths = [len(umi) for umi in umi_sequences]
        length_stats = {
            'mean': np.mean(umi_lengths),
            'std': np.std(umi_lengths),
            'min': min(umi_lengths),
            'max': max(umi_lengths),
            'median': np.median(umi_lengths)
        }
        
        # UMI composition analysis
        all_bases = ''.join(umi_sequences)
        base_counts = Counter(all_bases)
        total_bases = len(all_bases)
        composition = {base: count/total_bases for base, count in base_counts.items()}
        
        # UMI diversity (Shannon entropy)
        umi_counts = Counter(umi_sequences)
        total_unique = len(umi_counts)
        if total_unique > 1:
            proportions = [count/total_umis for count in umi_counts.values()]
            diversity = -sum(p * np.log2(p) for p in proportions if p > 0)
        else:
            diversity = 0
        
        # Read count analysis
        read_counts = {}
        if 'read_count' in df.columns:
            read_counts = {
                'mean': df['read_count'].mean(),
                'std': df['read_count'].std(),
                'min': df['read_count'].min(),
                'max': df['read_count'].max(),
                'median': df['read_count'].median()
            }
        
        return {
            'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
            'total_umis': total_umis,
            'unique_umis': unique_umis,
            'umi_length_stats': length_stats,
            'umi_composition': composition,
            'umi_diversity': diversity,
            'read_counts': read_counts
        }
        
    except Exception as e:
        logger.error(f"Error analyzing {umi_file}: {e}")
        return {
            'sample_name': umi_file.stem.replace('_Paired_UMIs_read_names', ''),
            'total_umis': 0,
            'unique_umis': 0,
            'umi_length_stats': {},
            'umi_composition': {},
            'umi_diversity': 0,
            'read_counts': {},
            'error': str(e)
        }


def umi_meta_analysis(umi_dir: str, output_dir: str, samples: Optional[List[str]] = None):
    """
    Perform meta-analysis on UMI mapping data
    
    Args:
        umi_dir: Directory containing UMI mapping CSV files
        output_dir: Output directory for results
        samples: Optional list of sample names
    """
    logger.info("Starting UMI meta-analysis")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find UMI files
    umi_path = Path(umi_dir)
    umi_files = list(umi_path.glob("*_Paired_UMIs_read_names.csv"))
    
    if not umi_files:
        raise ValueError(f"No UMI mapping files found in {umi_dir}")
    
    # Filter by samples if provided
    if samples:
        umi_files = [f for f in umi_files if any(sample in f.name for sample in samples)]
    
    logger.info(f"Found {len(umi_files)} UMI files to analyze")
    
    # Analyze each sample
    sample_results = []
    for umi_file in umi_files:
        result = analyze_umi_patterns(umi_file)
        sample_results.append(result)
    
    # Generate summary
    generate_umi_summary(sample_results, output_path)
    
    logger.info(f"UMI meta-analysis complete. Results saved to {output_dir}")


def generate_umi_summary(sample_results: List[Dict], output_path: Path):
    """
    Generate summary files from UMI analysis results
    
    Args:
        sample_results: List of UMI analysis results
        output_path: Path to output directory
    """
    logger.info("Generating UMI summary")
    
    # Create summary DataFrame
    summary_data = []
    for result in sample_results:
        summary_data.append({
            'Sample': result['sample_name'],
            'Total_UMIs': result['total_umis'],
            'Unique_UMIs': result['unique_umis'],
            'UMI_Diversity': result['umi_diversity'],
            'Mean_Length': result['umi_length_stats'].get('mean', 0),
            'Length_Std': result['umi_length_stats'].get('std', 0),
            'Min_Length': result['umi_length_stats'].get('min', 0),
            'Max_Length': result['umi_length_stats'].get('max', 0),
            'A_Composition': result['umi_composition'].get('A', 0),
            'T_Composition': result['umi_composition'].get('T', 0),
            'G_Composition': result['umi_composition'].get('G', 0),
            'C_Composition': result['umi_composition'].get('C', 0),
            'N_Composition': result['umi_composition'].get('N', 0)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_path / "umi_summary.csv", index=False)
    
    # Create detailed report
    report = {
        "analysis_metadata": {
            "total_samples": len(sample_results),
            "total_umis": sum(r['total_umis'] for r in sample_results),
            "total_unique_umis": sum(r['unique_umis'] for r in sample_results),
            "overall_diversity": np.mean([r['umi_diversity'] for r in sample_results if r['umi_diversity'] > 0])
        },
        "results": sample_results
    }
    
    with open(output_path / "umi_analysis_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info("UMI summary files generated successfully")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze UMI patterns and statistics")
    parser.add_argument("--umi-dir", required=True, help="Directory containing UMI mapping CSV files")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--samples", help="Comma-separated list of sample names")
    
    args = parser.parse_args()
    
    samples = None
    if args.samples:
        samples = [s.strip() for s in args.samples.split(',')]
    
    umi_meta_analysis(
        umi_dir=args.umi_dir,
        output_dir=args.output_dir,
        samples=samples
    )
