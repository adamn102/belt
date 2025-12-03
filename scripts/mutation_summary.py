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
from Bio import SeqIO

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


def load_reference_sequences(reference_dir: Path) -> Dict[str, int]:
    """
    Load reference sequences and return their lengths
    
    Args:
        reference_dir: Directory containing reference FASTA files
        
    Returns:
        Dictionary mapping reference base name to sequence length
    """
    reference_lengths = {}
    ref_path = Path(reference_dir)
    
    if not ref_path.exists():
        logger.warning(f"Reference directory not found: {reference_dir}")
        return reference_lengths
    
    # Find all FASTA files (exclude temp files and context CSVs)
    ref_files = [f for f in ref_path.glob("*.fa") 
                 if not f.name.startswith("~$") and not f.name.startswith(".")]
    ref_files.extend([f for f in ref_path.glob("*.fasta") 
                      if not f.name.startswith("~$") and not f.name.startswith(".")])
    
    for ref_file in ref_files:
        try:
            # Extract base name from filename (e.g., "LP249.fa" -> "LP249")
            base_name = ref_file.stem
            
            # Skip if this is a context file or other non-reference file
            if "_nucleotide_context" in base_name or base_name.endswith("_context"):
                continue
            
            # Read the first sequence from the FASTA file
            with open(ref_file, "r") as f:
                record = SeqIO.read(f, "fasta")
                seq_length = len(record.seq)
                reference_lengths[base_name] = seq_length
                logger.info(f"Loaded reference {base_name}: {seq_length} bases")
        except Exception as e:
            logger.error(f"Error loading reference from {ref_file}: {e}")
            continue
    
    return reference_lengths


def validate_mutation_analysis(sample_results: List[Dict], reference_dir: Optional[Path] = None):
    """
    Validate mutation analysis by asserting reference sequence lengths match expectations
    
    Args:
        sample_results: List of mutation analysis results from individual samples
        reference_dir: Optional directory containing reference FASTA files
        
    Raises:
        AssertionError: If validation fails
    """
    if not reference_dir or not reference_dir.exists():
        logger.warning(f"Reference directory not provided or doesn't exist, skipping validation")
        return
    
    logger.info("Validating mutation analysis against reference sequences")
    
    # Load reference sequence lengths
    reference_lengths = load_reference_sequences(reference_dir)
    
    if not reference_lengths:
        logger.warning("No reference sequences found, skipping validation")
        return
    
    # Validate each sample
    for result in sample_results:
        sample_name = result['sample_name']
        reference_base = sample_name.split('-')[0]
        
        if reference_base not in reference_lengths:
            logger.warning(f"Reference {reference_base} not found for sample {sample_name}, skipping validation")
            continue
        
        expected_length = reference_lengths[reference_base]
        
        # Check that mutations don't exceed reference length
        if 'mutations_by_position' in result:
            max_position = -1
            for mutation in result['mutations_by_position']:
                position = mutation.get('position', -1)
                if position > max_position:
                    max_position = position
            
            if max_position >= 0:
                # Reference length is 0-indexed, so max valid position is length - 1
                max_valid_position = expected_length - 1
                assert max_position <= max_valid_position, (
                    f"Sample {sample_name}: Mutation at position {max_position} exceeds "
                    f"reference length {expected_length} (max valid position: {max_valid_position})"
                )
                logger.debug(f"Sample {sample_name}: All mutations within reference bounds "
                           f"(max mutation position: {max_position}, reference length: {expected_length})")
        
        # Log reference length for verification
        logger.info(f"Sample {sample_name}: Reference {reference_base} length = {expected_length} bases")
        
        # Assert that if we have nucleotide context, it matches the reference length
        # (This validates that the context file was generated from the correct reference)
        # Note: This check is done implicitly through the position bounds check above,
        # but we log it for clarity
        logger.debug(f"Sample {sample_name}: Reference length validation passed "
                    f"(reference: {expected_length} bases, all mutations within bounds)")
    
    logger.info("Mutation analysis validation complete")


def load_nucleotide_contexts(context_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load nucleotide context CSV files for all references
    
    Args:
        context_dir: Directory containing nucleotide context CSV files
        
    Returns:
        Dictionary mapping reference base name to context DataFrame
    """
    context_dict = {}
    context_path = Path(context_dir)
    
    if not context_path.exists():
        logger.warning(f"Nucleotide context directory not found: {context_dir}")
        return context_dict
    
    # Find all nucleotide context CSV files (exclude temp files)
    context_files = [f for f in context_path.glob("*_nucleotide_context.csv") 
                     if not f.name.startswith("~$") and not f.name.startswith(".")]
    
    if not context_files:
        logger.warning(f"No nucleotide context files found in {context_dir}")
        return context_dict
    
    for context_file in context_files:
        try:
            # Extract reference base name from filename (e.g., "LP249_nucleotide_context.csv" -> "LP249")
            base_name = context_file.stem.replace("_nucleotide_context", "")
            context_df = pd.read_csv(context_file)
            # Set position as index for faster lookup
            context_df.set_index('position', inplace=True)
            context_dict[base_name] = context_df
            logger.info(f"Loaded nucleotide context for {base_name}: {len(context_df)} positions")
        except Exception as e:
            logger.error(f"Error loading nucleotide context from {context_file}: {e}")
            continue
    
    return context_dict


def add_nucleotide_context_to_mutations(mutations: List[Dict], reference_base: str, 
                                         context_dict: Dict[str, pd.DataFrame]) -> List[Dict]:
    """
    Add nucleotide context information to mutations
    
    Args:
        mutations: List of mutation dictionaries
        reference_base: Base name of the reference (e.g., "LP249")
        context_dict: Dictionary mapping reference base to context DataFrame
        
    Returns:
        List of mutations with added context information
    """
    if reference_base not in context_dict:
        logger.warning(f"No nucleotide context found for reference {reference_base}")
        return mutations
    
    context_df = context_dict[reference_base]
    enriched_mutations = []
    
    for mutation in mutations:
        enriched_mutation = mutation.copy()
        position = mutation.get('position')
        
        # Handle position (should be 0-indexed, but check if it's in the context)
        if position in context_df.index:
            context_row = context_df.loc[position]
            enriched_mutation['nucleotide_context'] = {
                'tri_nucleotide_context': context_row['tri_nucleotide_context'],
                'nucleotide': context_row['nucleotide'],
                'upstream_nucleotide': context_row['upstream_nucleotide'],
                'downstream_nucleotide': context_row['downstream_nucleotide']
            }
        else:
            # Position not found in context (might be out of range or indexing issue)
            logger.debug(f"Position {position} not found in context for {reference_base}")
            enriched_mutation['nucleotide_context'] = None
        
        enriched_mutations.append(enriched_mutation)
    
    return enriched_mutations


def mutation_summary(sample_results: List[Dict], output_dir: Path, 
                     context_dir: Optional[Path] = None,
                     reference_dir: Optional[Path] = None):
    """
    Generate summary analysis from individual sample results
    
    Args:
        sample_results: List of mutation analysis results from individual samples
        output_dir: Output directory for summary files
        context_dir: Optional directory containing nucleotide context CSV files
        reference_dir: Optional directory containing reference FASTA files for validation
    """
    logger.info("Generating mutation summary analysis")
    
    # Validate mutation analysis against reference sequences
    if reference_dir:
        validate_mutation_analysis(sample_results, Path(reference_dir))
    
    # Load nucleotide contexts if provided
    context_dict = {}
    if context_dir:
        context_dict = load_nucleotide_contexts(Path(context_dir))
        
        # Validate that context files match reference lengths if reference_dir is provided
        if reference_dir:
            reference_lengths = load_reference_sequences(Path(reference_dir))
            for ref_base, context_df in context_dict.items():
                if ref_base in reference_lengths:
                    expected_length = reference_lengths[ref_base]
                    actual_positions = len(context_df)
                    assert actual_positions == expected_length, (
                        f"Nucleotide context for {ref_base} has {actual_positions} positions, "
                        f"but reference sequence has {expected_length} bases"
                    )
                    logger.info(f"Validated nucleotide context for {ref_base}: "
                              f"{actual_positions} positions match reference length")
    
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
    
    # Enrich mutations with nucleotide context if available
    enriched_results = []
    for result in sample_results:
        enriched_result = result.copy()
        sample_name = result['sample_name']
        
        # Extract reference base name from sample name (e.g., "LP249-NL43-WT-1-A" -> "LP249")
        reference_base = sample_name.split('-')[0]
        
        # Add context to mutations_by_position
        if 'mutations_by_position' in result and context_dict:
            enriched_result['mutations_by_position'] = add_nucleotide_context_to_mutations(
                result['mutations_by_position'],
                reference_base,
                context_dict
            )
            logger.info(f"Added nucleotide context to {len(enriched_result['mutations_by_position'])} mutations for {sample_name}")
        
        enriched_results.append(enriched_result)
    
    # Create detailed report
    report = {
        "analysis_metadata": {
            "total_samples": len(enriched_results),
            "total_consensus_sequences": sum(r['consensus_count'] for r in enriched_results),
            "total_mutations": sum(r['total_mutations'] for r in enriched_results),
            "total_length": sum(r['total_length'] for r in enriched_results),
            "overall_mutation_rate": sum(r['total_mutations'] for r in enriched_results) / sum(r['total_length'] for r in enriched_results) if sum(r['total_length'] for r in enriched_results) > 0 else 0,
            "nucleotide_context_enabled": len(context_dict) > 0
        },
        "results": enriched_results
    }
    
    with open(output_dir / "mutation_analysis_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    # Generate mutation distribution analysis
    all_mutations = []
    for result in enriched_results:
        for mutation in result['mutations_by_position']:
            # Add sample_id if not already present (for backward compatibility)
            if 'sample_id' not in mutation:
                mutation['sample_id'] = result['sample_name']
            all_mutations.append(mutation)
    
    if all_mutations:
        mutation_df = pd.DataFrame(all_mutations)
        mutation_df.to_csv(output_dir / "all_mutations.csv", index=False)
        
        # Position-based analysis
        position_stats = mutation_df.groupby('position').agg({
            'type': 'count',
            'sample_id': 'nunique'
        }).rename(columns={'type': 'mutation_count', 'sample_id': 'sample_count'})
        position_stats.to_csv(output_dir / "mutation_position_stats.csv")
    else:
        # Create empty files when no mutations are found
        pd.DataFrame(columns=['position', 'type', 'ref_base', 'cons_base', 'umi_id', 'sample_id']).to_csv(output_dir / "all_mutations.csv", index=False)
        pd.DataFrame(columns=['mutation_count', 'sample_count']).to_csv(output_dir / "mutation_position_stats.csv", index=True)
    
    logger.info("Summary analysis complete")


def create_mutation_summary_from_files(output_dir: str, input_dir: Optional[str] = None,
                                       context_dir: Optional[str] = None,
                                       reference_dir: Optional[str] = None):
    """
    Create summary analysis from individual sample mutation files
    
    Args:
        output_dir: Output directory for summary results
        input_dir: Input directory containing sample results (optional)
        context_dir: Directory containing nucleotide context CSV files (optional)
        reference_dir: Directory containing reference FASTA files for validation (optional)
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
    
    # Generate summary with optional context and reference validation
    context_path = Path(context_dir) if context_dir else None
    reference_path = Path(reference_dir) if reference_dir else None
    mutation_summary(sample_results, output_path, context_dir=context_path, reference_dir=reference_path)
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
    parser.add_argument("--context-dir", help="Directory containing nucleotide context CSV files")
    parser.add_argument("--reference-dir", help="Directory containing reference FASTA files for validation")
    
    args = parser.parse_args()
    
    create_mutation_summary_from_files(
        output_dir=args.output_dir,
        input_dir=args.input_dir,
        context_dir=args.context_dir,
        reference_dir=args.reference_dir
    )
