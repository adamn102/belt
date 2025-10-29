#!/usr/bin/env python3
"""
BELT Mutation Summary Module
Generates summary analysis from individual sample mutation results.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional
import glob
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def generate_consensus_sequence_data(sample_results: List[Dict], output_dir: Path):
    """
    Generate detailed consensus sequence data CSV
    
    Args:
        sample_results: List of mutation analysis results from individual samples
        output_dir: Output directory for summary files
    """
    logger.info("Generating consensus sequence data CSV")
    
    consensus_data = []
    for result in sample_results:
        sample_name = result['sample_name']
        
        # Process each consensus sequence in the sample
        for read_analysis in result.get('per_read_analysis', []):
            consensus_id = read_analysis['consensus_id']
            sequence_length = read_analysis['sequence_length']
            mutation_count = read_analysis['mutation_count']
            mutation_rate = read_analysis.get('mutation_rate', 0)
            
            # Count different mutation types for this consensus
            substitutions = 0
            insertions = 0
            deletions = 0
            
            for mutation in read_analysis.get('mutations', []):
                if mutation['type'] == 'substitution':
                    substitutions += 1
                elif mutation['type'] == 'insertion':
                    insertions += 1
                elif mutation['type'] == 'deletion':
                    deletions += 1
            
            consensus_data.append({
                'Sample': sample_name,
                'UMI_ID': consensus_id,
                'Total_Bases': sequence_length,
                'Total_Mutations': mutation_count,
                'Mutation_Rate': mutation_rate,
                'Substitutions': substitutions,
                'Insertions': insertions,
                'Deletions': deletions
            })
    
    # Create DataFrame and save
    if consensus_data:
        consensus_df = pd.DataFrame(consensus_data)
        consensus_df.to_csv(output_dir / "consensus_sequence_data.csv", index=False)
        logger.info(f"Saved consensus sequence data for {len(consensus_data)} sequences to consensus_sequence_data.csv")
    else:
        # Create empty file if no data
        empty_df = pd.DataFrame(columns=[
            'Sample', 'UMI_ID', 'Total_Bases', 'Total_Mutations', 
            'Mutation_Rate', 'Substitutions', 'Insertions', 'Deletions'
        ])
        empty_df.to_csv(output_dir / "consensus_sequence_data.csv", index=False)
        logger.info("Created empty consensus_sequence_data.csv (no consensus sequences found)")


def mutation_summary(sample_results: List[Dict], output_dir: Path):
    """
    Generate summary analysis from individual sample results
    
    Args:
        sample_results: List of mutation analysis results from individual samples
        output_dir: Output directory for summary files
    """
    logger.info("Generating mutation summary analysis")
    
    # Create summary DataFrame
    summary_data = []
    for result in sample_results:
        summary_data.append({
            'Sample': result['sample_name'],
            'Consensus_Count': result['consensus_count'],
            'Total_Length': result['total_length'],
            'Total_Mutations': result['total_mutations'],
            'Average_Mutation_Rate': result['average_mutation_rate'],
            'Substitutions': result['mutation_types']['substitution'],
            'Insertions': result['mutation_types']['insertion'],
            'Deletions': result['mutation_types']['deletion']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / "mutation_summary.csv", index=False)
    
    # Generate consensus sequence data CSV
    generate_consensus_sequence_data(sample_results, output_dir)
    
    # Create detailed report
    report = {
        "analysis_metadata": {
            "total_samples": len(sample_results),
            "total_consensus_sequences": sum(r['consensus_count'] for r in sample_results),
            "total_mutations": sum(r['total_mutations'] for r in sample_results),
            "total_length": sum(r['total_length'] for r in sample_results),
            "overall_mutation_rate": sum(r['total_mutations'] for r in sample_results) / sum(r['total_length'] for r in sample_results) if sum(r['total_length'] for r in sample_results) > 0 else 0
        },
        "results": sample_results
    }
    
    with open(output_dir / "mutation_analysis_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    # Generate mutation distribution analysis
    all_mutations = []
    for result in sample_results:
        for mutation in result['mutations_by_position']:
            mutation['sample'] = result['sample_name']
            all_mutations.append(mutation)
    
    if all_mutations:
        mutation_df = pd.DataFrame(all_mutations)
        mutation_df.to_csv(output_dir / "all_mutations.csv", index=False)
        
        # Position-based analysis
        position_stats = mutation_df.groupby('position').agg({
            'type': 'count',
            'sample': 'nunique'
        }).rename(columns={'type': 'mutation_count', 'sample': 'sample_count'})
        position_stats.to_csv(output_dir / "mutation_position_stats.csv")
    else:
        # Create empty files when no mutations are found
        pd.DataFrame(columns=['position', 'type', 'ref_base', 'cons_base', 'umi_id', 'sample']).to_csv(output_dir / "all_mutations.csv", index=False)
        pd.DataFrame(columns=['mutation_count', 'sample_count']).to_csv(output_dir / "mutation_position_stats.csv", index=True)
    
    logger.info("Summary analysis complete")


def create_mutation_summary_from_files(output_dir: str, input_dir: Optional[str] = None):
    """
    Create summary analysis from individual sample mutation files
    
    Args:
        output_dir: Output directory for summary results
        input_dir: Input directory containing sample results (optional)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Use input_dir if provided, otherwise look for step 9 results
    if input_dir:
        search_dir = input_dir
    else:
        # Look for step 9 results in the parent directory
        search_dir = str(output_path.parent / "9-mutation_calling")
    
    # Find all sample result files
    sample_files = glob.glob(f"{search_dir}/*_mutations.json")
    
    if not sample_files:
        logger.error(f"No sample result files found in {search_dir}")
        logger.error(f"Expected files like: {search_dir}/*_mutations.json")
        return
    
    logger.info(f"Found {len(sample_files)} sample result files")
    
    # Load all sample results
    sample_results = []
    for sample_file in sample_files:
        try:
            with open(sample_file, 'r') as f:
                sample_results.append(json.load(f))
        except Exception as e:
            logger.error(f"Error loading {sample_file}: {e}")
            continue
    
    if not sample_results:
        logger.error("No valid sample results loaded")
        return
    
    # Generate summary
    mutation_summary(sample_results, output_path)
    logger.info("Mutation summary analysis complete")


def analyze_mutations(consensus_dir: str, reference_dir: str, output_dir: str, 
                     samples: Optional[List[str]] = None, threads: int = 4):
    """
    Main function to analyze mutations for all samples (legacy compatibility)
    
    Args:
        consensus_dir: Directory containing consensus FASTA files
        reference_dir: Directory containing reference FASTA files
        output_dir: Output directory for results
        samples: Optional list of sample names
        threads: Number of threads for parallel processing
    """
    logger.info("Starting mutation analysis for all samples")
    
    # Import mutation_caller here to avoid circular imports
    from mutation_caller import mutation_caller
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find consensus files
    consensus_path = Path(consensus_dir)
    consensus_files = list(consensus_path.glob("*-compiled-consensus.fa"))
    
    if not consensus_files:
        raise ValueError(f"No consensus files found in {consensus_dir}")
    
    # Find reference files
    reference_path = Path(reference_dir)
    reference_files = list(reference_path.glob("*.fa"))
    
    if not reference_files:
        raise ValueError(f"No reference files found in {reference_dir}")
    
    # Process each consensus file
    sample_results = []
    for consensus_file in consensus_files:
        sample_name = consensus_file.stem.replace("-compiled-consensus", "")
        
        # Find matching reference file
        reference_file = None
        for ref_file in reference_files:
            if sample_name in ref_file.name or "LP249" in ref_file.name:
                reference_file = ref_file
                break
        
        if not reference_file:
            logger.warning(f"No reference file found for sample {sample_name}")
            continue
        
        logger.info(f"Found reference file for {sample_name}: {reference_file}")
        
        try:
            result = mutation_caller(consensus_file, reference_file, sample_name, output_path)
            sample_results.append(result)
        except Exception as e:
            logger.error(f"Error analyzing {sample_name}: {e}")
            continue
    
    # Generate summary
    mutation_summary(sample_results, output_path)
    
    logger.info(f"Mutation analysis complete. Results saved to {output_dir}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate mutation summary analysis")
    parser.add_argument("--output-dir", required=True, help="Output directory for summary results")
    parser.add_argument("--input-dir", help="Input directory containing sample results (default: data/9-mutation_calling)")
    
    args = parser.parse_args()
    
    create_mutation_summary_from_files(
        output_dir=args.output_dir,
        input_dir=args.input_dir
    )
