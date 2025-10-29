#!/usr/bin/env python3
"""
BELT Mutation Caller Module
Analyzes mutations for individual samples by comparing consensus sequences to reference sequences.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging

# BioPython imports for sequence analysis
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_sequences(fasta_path: Path) -> Dict[str, SeqRecord]:
    """Load sequences from a FASTA file and return as dictionary"""
    sequences = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences[record.id] = record
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
        return sequences
    except Exception as e:
        logger.error(f"Error loading sequences from {fasta_path}: {e}")
        return {}


def align_sequences(consensus_seq: Seq, reference_seq: Seq) -> Tuple[str, str, Dict]:
    """Align consensus sequence to reference and return alignment details"""
    try:
        # Create aligner with appropriate parameters for DNA
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Global alignment
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        
        # Perform alignment
        alignments = aligner.align(reference_seq, consensus_seq)
        
        if not alignments:
            return "", "", {}
        
        # Get the best alignment
        best_alignment = alignments[0]
        
        # Extract aligned sequences
        ref_aligned = str(best_alignment[0])
        cons_aligned = str(best_alignment[1])
        
        # Calculate alignment statistics
        alignment_stats = {
            'score': best_alignment.score,
            'length': len(ref_aligned),
            'matches': sum(1 for r, c in zip(ref_aligned, cons_aligned) if r == c and r != '-'),
            'mismatches': sum(1 for r, c in zip(ref_aligned, cons_aligned) if r != c and r != '-' and c != '-'),
            'gaps': sum(1 for r, c in zip(ref_aligned, cons_aligned) if r == '-' or c == '-')
        }
        
        return ref_aligned, cons_aligned, alignment_stats
        
    except Exception as e:
        logger.error(f"Error aligning sequences: {e}")
        return "", "", {}


def find_mutations(ref_aligned: str, cons_aligned: str) -> List[Dict]:
    """
    Find mutations between aligned sequences.
    Groups consecutive indels together as single events with size information.
    """
    mutations = []
    ref_pos = 0
    i = 0
    
    while i < len(ref_aligned):
        ref_char = ref_aligned[i]
        cons_char = cons_aligned[i]
        
        if ref_char == '-' and cons_char != '-':
            # Start of insertion - collect consecutive insertions
            insertion_start = i
            insertion_bases = []
            insertion_ref_pos = ref_pos
            
            while i < len(ref_aligned) and ref_aligned[i] == '-' and cons_aligned[i] != '-':
                insertion_bases.append(cons_aligned[i])
                i += 1
            
            # Record as single insertion with size
            mutations.append({
                "type": "insertion",
                "ref_position": insertion_ref_pos,
                "ref_base": "-",
                "cons_base": "".join(insertion_bases) if len(insertion_bases) > 1 else insertion_bases[0],
                "size": len(insertion_bases),
                "inserted_sequence": "".join(insertion_bases)
            })
            continue
            
        elif ref_char != '-' and cons_char == '-':
            # Start of deletion - collect consecutive deletions
            deletion_start_pos = ref_pos
            deletion_bases = []
            
            while i < len(ref_aligned) and ref_aligned[i] != '-' and cons_aligned[i] == '-':
                deletion_bases.append(ref_aligned[i])
                ref_pos += 1
                i += 1
            
            # Record as single deletion with size
            mutations.append({
                "type": "deletion",
                "ref_position": deletion_start_pos,
                "ref_base": "".join(deletion_bases) if len(deletion_bases) > 1 else deletion_bases[0],
                "cons_base": "-",
                "size": len(deletion_bases),
                "deleted_sequence": "".join(deletion_bases)
            })
            continue
            
        elif ref_char != '-' and cons_char != '-' and ref_char != cons_char:
            # Substitution
            mutations.append({
                "type": "substitution",
                "ref_position": ref_pos,
                "ref_base": ref_char,
                "cons_base": cons_char
            })
            ref_pos += 1
            
        elif ref_char != '-' and cons_char != '-':
            # Match
            ref_pos += 1
        
        i += 1
    
    return mutations


def mutation_caller(consensus_file: Path, reference_file: Path, sample_name: str, 
                   output_dir: Path) -> Dict:
    """
    Analyze mutations for a single sample and save results
    
    Args:
        consensus_file: Path to consensus FASTA file
        reference_file: Path to reference FASTA file
        sample_name: Name of the sample
        output_dir: Output directory for individual sample results
    
    Returns:
        Dictionary containing mutation analysis results for the sample
    """
    logger.info(f"Analyzing mutations for sample: {sample_name}")
    
    # Load consensus sequences
    consensus_sequences = load_sequences(consensus_file)
    if not consensus_sequences:
        raise ValueError(f"No consensus sequences found in {consensus_file}")
    
    # Load reference sequence
    reference_sequences = load_sequences(reference_file)
    if not reference_sequences:
        raise ValueError(f"No reference sequences found in {reference_file}")
    
    # Use the first reference sequence
    reference_id = list(reference_sequences.keys())[0]
    reference_seq = reference_sequences[reference_id].seq
    
    # Analyze mutations for each consensus sequence
    mutations_by_position = []
    mutation_types = {"substitution": 0, "insertion": 0, "deletion": 0}
    total_mutations = 0
    total_length = 0
    per_read_analysis = []
    
    for cons_id, cons_record in consensus_sequences.items():
        cons_seq = cons_record.seq
        
        # Perform alignment
        ref_aligned, cons_aligned, alignment_stats = align_sequences(cons_seq, reference_seq)
        
        if not ref_aligned or not cons_aligned:
            logger.warning(f"Failed to align sequence {cons_id}")
            continue
        
        # Find mutations
        mutations = find_mutations(ref_aligned, cons_aligned)
        
        # Process mutations
        read_mutations = []
        for mutation in mutations:
            mutation_types[mutation["type"]] += 1
            mutations_by_position.append({
                "position": mutation["ref_position"],
                "type": mutation["type"],
                "ref_base": mutation["ref_base"],
                "cons_base": mutation["cons_base"],
                "umi_id": cons_id
            })
            read_mutations.append(mutation)
        
        # Per-read analysis
        per_read_analysis.append({
            "consensus_id": cons_id,
            "sequence_length": len(cons_seq),
            "mutation_count": len(mutations),
            "mutation_rate": len(mutations) / len(cons_seq) if len(cons_seq) > 0 else 0,
            "mutations": read_mutations
        })
        
        total_mutations += len(mutations)
        total_length += len(cons_seq)
    
    # Calculate statistics
    consensus_count = len(consensus_sequences)
    average_mutation_rate = total_mutations / total_length if total_length > 0 else 0
    
    # Create sample result
    sample_result = {
        "sample_name": sample_name,
        "consensus_count": consensus_count,
        "total_length": total_length,
        "total_mutations": total_mutations,
        "average_mutation_rate": average_mutation_rate,
        "mutation_types": mutation_types,
        "mutations_by_position": mutations_by_position,
        "per_read_analysis": per_read_analysis
    }
    
    # Save individual sample results
    sample_output_dir = output_dir / "samples" / sample_name
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save detailed results
    with open(sample_output_dir / f"{sample_name}_mutations.json", 'w') as f:
        json.dump(sample_result, f, indent=2)
    
    # Save per-read analysis as CSV
    per_read_df = pd.DataFrame([
        {
            'Consensus_ID': read['consensus_id'],
            'Sequence_Length': read['sequence_length'],
            'Mutation_Count': read['mutation_count'],
            'Mutation_Rate': read['mutation_rate']
        }
        for read in per_read_analysis
    ])
    per_read_df.to_csv(sample_output_dir / f"{sample_name}_per_read.csv", index=False)
    
    logger.info(f"Sample {sample_name} analysis complete: {consensus_count} sequences, {total_mutations} mutations")
    
    return sample_result


def call_mutations_for_sample(consensus_file: str, reference_file: str, 
                              sample_name: str, output_dir: str):
    """
    Wrapper function for calling mutations for a single sample (CLI interface)
    
    Args:
        consensus_file: Path to consensus FASTA file
        reference_file: Path to reference FASTA file
        sample_name: Name of the sample
        output_dir: Output directory for results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    result = mutation_caller(
        consensus_file=Path(consensus_file),
        reference_file=Path(reference_file),
        sample_name=sample_name,
        output_dir=output_path
    )
    
    # Save the result in the expected format for step 9
    output_file = output_path / f"{sample_name}_mutations.json"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"Completed mutation analysis for {sample_name}")
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Call mutations for a single sample")
    parser.add_argument("--consensus-file", required=True, help="Path to consensus FASTA file")
    parser.add_argument("--reference-file", required=True, help="Path to reference FASTA file")
    parser.add_argument("--sample-name", required=True, help="Name of the sample")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    
    args = parser.parse_args()
    
    call_mutations_for_sample(
        consensus_file=args.consensus_file,
        reference_file=args.reference_file,
        sample_name=args.sample_name,
        output_dir=args.output_dir
    )
