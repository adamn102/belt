"""
Improved consensus calling with pre-filtering for large deletions.

This module provides alternative consensus calling methods that handle
size-variable large deletions better than standard bcftools.
"""

import pysam
import subprocess
import shlex
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import shutil


def detect_deletion_regions(bam_path, umi, min_deletion_size=50, min_read_support=2):
    """
    Detect large deletion regions from CIGAR strings.
    
    Returns list of (position, sizes_dict) where sizes_dict maps
    deletion_size -> read_count
    """
    deletion_regions = {}
    
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if not read.has_tag('ZU'):
                continue
            if read.get_tag('ZU') != umi:
                continue
            
            ref_pos = read.reference_start
            for op, length in read.cigar:
                if op == 2 and length >= min_deletion_size:  # Large deletion
                    if ref_pos not in deletion_regions:
                        deletion_regions[ref_pos] = Counter()
                    deletion_regions[ref_pos][length] += 1
                
                if op in [0, 2, 3, 7, 8]:  # Consume reference
                    ref_pos += length
    
    # Filter to regions with sufficient support
    filtered = {
        pos: sizes for pos, sizes in deletion_regions.items()
        if sum(sizes.values()) >= min_read_support
    }
    
    return filtered


def filter_reads_by_deletion_size(bam_path, umi, deletion_pos, size_tolerance=50, output_bam=None):
    """
    Filter reads to those with deletion sizes within tolerance of most common size.
    
    Returns: (filtered_bam_path, most_common_size, read_count)
    """
    # First pass: find deletion sizes at target position
    deletion_sizes = []
    read_deletion_map = {}  # read_name -> deletion_size
    
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        header = bam.header
        for read in bam.fetch(until_eof=True):
            if not read.has_tag('ZU') or read.get_tag('ZU') != umi:
                continue
            
            ref_pos = read.reference_start
            read_del_size = None
            
            for op, length in read.cigar:
                if op == 2:  # Deletion
                    # Check if deletion is near target position
                    if abs(ref_pos - deletion_pos) < 10:
                        deletion_sizes.append(length)
                        read_del_size = length
                if op in [0, 2, 3, 7, 8]:
                    ref_pos += length
            
            if read_del_size is not None:
                read_deletion_map[read.query_name] = read_del_size
    
    if not deletion_sizes:
        return None, None, 0
    
    # Find most common deletion size
    size_counts = Counter(deletion_sizes)
    most_common_size = size_counts.most_common(1)[0][0]
    
    # Filter reads within tolerance
    if output_bam is None:
        output_bam = str(Path(bam_path).parent / f"{umi}_filtered.bam")
    
    filtered_count = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam_in, \
         pysam.AlignmentFile(output_bam, 'wb', header=bam_in.header) as bam_out:
        
        for read in bam_in.fetch(until_eof=True):
            if not read.has_tag('ZU') or read.get_tag('ZU') != umi:
                continue
            
            read_del_size = read_deletion_map.get(read.query_name)
            if read_del_size is not None:
                if abs(read_del_size - most_common_size) <= size_tolerance:
                    bam_out.write(read)
                    filtered_count += 1
    
    # Index the filtered BAM
    pysam.index(output_bam)
    
    return output_bam, most_common_size, filtered_count


def umi_consensus_with_prefiltering(umi, bam, ref, threads=1, 
                                   deletion_size_tolerance=50,
                                   use_filtered=True, **kwargs):
    """
    Consensus calling with pre-filtering for size-variable deletions.
    
    Strategy:
    1. Detect large deletion regions from CIGAR strings
    2. If size variation detected, filter reads to most common size
    3. Run standard consensus calling on filtered reads
    4. Fall back to unfiltered if filtering removes too many reads
    """
    from .consensus import umi_consensus_bcftools
    
    bam_path = Path(bam)
    ref_path = Path(ref)
    
    # Detect deletion regions
    deletion_regions = detect_deletion_regions(str(bam_path), umi, min_deletion_size=50)
    
    if not deletion_regions or not use_filtered:
        # No large deletions detected or filtering disabled, use standard method
        return umi_consensus_bcftools(umi, bam, ref, threads=threads, **kwargs)
    
    # Find region with most size variation
    max_variation = 0
    target_pos = None
    target_sizes = None
    
    for pos, sizes in deletion_regions.items():
        if len(sizes) > 1:  # Multiple sizes
            size_range = max(sizes.keys()) - min(sizes.keys())
            if size_range > max_variation:
                max_variation = size_range
                target_pos = pos
                target_sizes = sizes
    
    if target_pos is None or max_variation < 50:
        # No significant size variation, use standard method
        return umi_consensus_bcftools(umi, bam, ref, threads=threads, **kwargs)
    
    # Filter reads
    print(f"Detected size-variable deletion at position {target_pos}")
    print(f"  Size range: {min(target_sizes.keys())}-{max(target_sizes.keys())}bp")
    print(f"  Filtering reads to most common size ±{deletion_size_tolerance}bp...")
    
    filtered_bam, most_common_size, filtered_count = filter_reads_by_deletion_size(
        str(bam_path), umi, target_pos, deletion_size_tolerance
    )
    
    if filtered_bam is None or filtered_count < 3:
        print(f"  Warning: Filtering resulted in <3 reads, using unfiltered BAM")
        return umi_consensus_bcftools(umi, bam, ref, threads=threads, **kwargs)
    
    print(f"  Filtered to {filtered_count} reads with deletion size {most_common_size}±{deletion_size_tolerance}bp")
    
    # Run consensus on filtered BAM
    try:
        result = umi_consensus_bcftools(
            umi, filtered_bam, ref, threads=threads, 
            skip_extraction=True, **kwargs
        )
        return result
    finally:
        # Clean up filtered BAM
        try:
            Path(filtered_bam).unlink(missing_ok=True)
            Path(f"{filtered_bam}.bai").unlink(missing_ok=True)
        except:
            pass


def manual_consensus_with_deletion(bam_path, umi, ref_path, deletion_pos, 
                                  consensus_del_size=None):
    """
    Manually construct consensus by taking majority deletion size from CIGAR strings.
    
    This bypasses variant calling entirely and directly constructs the consensus.
    """
    # Read reference
    with open(ref_path) as f:
        ref_record = SeqIO.read(f, "fasta")
        ref_seq = str(ref_record.seq)
    
    # Collect deletion sizes at target position
    deletion_sizes = []
    
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if not read.has_tag('ZU') or read.get_tag('ZU') != umi:
                continue
            
            ref_pos = read.reference_start
            for op, length in read.cigar:
                if op == 2:  # Deletion
                    if abs(ref_pos - deletion_pos) < 10:  # Near target
                        deletion_sizes.append(length)
                if op in [0, 2, 3, 7, 8]:
                    ref_pos += length
    
    if not deletion_sizes:
        # No deletion detected, return reference
        return ref_seq, 0
    
    # Use most common deletion size
    if consensus_del_size is None:
        size_counts = Counter(deletion_sizes)
        consensus_del_size = size_counts.most_common(1)[0][0]
    
    # Construct consensus: reference with deletion
    if deletion_pos + consensus_del_size <= len(ref_seq):
        consensus_seq = ref_seq[:deletion_pos] + ref_seq[deletion_pos + consensus_del_size:]
    else:
        consensus_seq = ref_seq[:deletion_pos]
    
    return consensus_seq, consensus_del_size


def hybrid_consensus_calling(umi, bam, ref, threads=1, **kwargs):
    """
    Hybrid approach: Use bcftools for most of the sequence, but manually
    handle known problematic deletion regions.
    """
    from .consensus import umi_consensus_bcftools
    
    # First, try standard consensus calling
    try:
        consensus_file = umi_consensus_bcftools(umi, bam, ref, threads=threads, **kwargs)
        
        # Check if consensus has expected length (indicates deletion was called)
        with open(consensus_file) as f:
            consensus_lines = f.readlines()
            consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
        
        with open(ref) as f:
            ref_record = SeqIO.read(f, "fasta")
            ref_seq = str(ref_record.seq)
        
        length_diff = len(ref_seq) - len(consensus_seq)
        
        # If consensus is significantly shorter, deletion was likely called
        if length_diff > 50:
            return consensus_file
        
        # Otherwise, try pre-filtering approach
        print("Standard consensus may have missed deletion, trying pre-filtering...")
        return umi_consensus_with_prefiltering(umi, bam, ref, threads=threads, **kwargs)
        
    except Exception as e:
        print(f"Standard consensus failed: {e}")
        print("Trying pre-filtering approach...")
        return umi_consensus_with_prefiltering(umi, bam, ref, threads=threads, **kwargs)

