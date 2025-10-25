#!/usr/bin/env python3
"""
Generate synthetic demo data for BELT pipeline
This script creates small, GitHub-friendly test data including:
- Reference sequence with UMI motif positions
- FASTQ reads with realistic UMI barcodes
- UMI whitelist

Run this script to regenerate demo data if needed.
"""

import random
from pathlib import Path

# Set seed for reproducibility
random.seed(42)

def generate_dna(length):
    """Generate random DNA sequence"""
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=length))

def add_mutation(seq, rate=0.02):
    """Add random mutations to sequence"""
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice(['A', 'C', 'G', 'T'])
    return ''.join(seq)

def generate_quality_scores(length, mean_qual=30):
    """Generate random quality scores"""
    quals = []
    for _ in range(length):
        q = max(10, min(40, int(random.gauss(mean_qual, 5))))
        quals.append(chr(q + 33))
    return ''.join(quals)

def main():
    print("Generating BELT demo data...")
    print("=" * 70)
    
    # Create reference with UMI motif markers
    # Structure: [prefix] NNNNNNNNNNNN [middle] NNNNNNNNNNNN [suffix]
    umi_motif = "N" * 12
    prefix_len = 500
    middle_len = 1400
    suffix_len = 500
    
    prefix_seq = generate_dna(prefix_len)
    middle_seq = generate_dna(middle_len)
    suffix_seq = generate_dna(suffix_len)
    
    # Build reference with N's marking UMI positions
    ref_seq = prefix_seq + umi_motif + middle_seq + umi_motif + suffix_seq
    
    ref_path = Path("data/0-reference_sequence/demo_ref.fa")
    ref_path.write_text(f">demo_reference\n{ref_seq}\n")
    print(f"✅ Reference created: {ref_path} ({len(ref_seq)} bp)")
    print(f"   UMI1 position: {prefix_len}")
    print(f"   UMI2 position: {prefix_len + 12 + middle_len}")
    
    # Also copy for consensus reference
    cons_ref_path = Path("data/0-consensus_reference/demo_ref.fa")
    cons_ref_path.write_text(f">demo_reference\n{ref_seq}\n")
    print(f"✅ Consensus reference created: {cons_ref_path}")
    
    # Generate whitelist with 12 UMIs
    umis = []
    for i in range(12):
        umi = ''.join(random.choices(['A', 'C', 'G', 'T'], k=12))
        umis.append(umi)
    
    whitelist_path = Path("data/0-library_whitelist/demo_whitelist.csv")
    whitelist_path.write_text("UMI1,UMI2\n" + "\n".join([f"{umi},{umi}" for umi in umis]))
    print(f"✅ Whitelist created: {whitelist_path} ({len(umis)} UMI pairs)")
    
    # Generate synthetic FASTQ with 150 reads
    reads = []
    umi_distribution = {umi: 0 for umi in umis}
    
    for read_id in range(150):
        # Select a random UMI pair (biased distribution for realistic data)
        if read_id < 80:  # First 80 reads use first 6 UMIs heavily
            umi1, umi2 = random.choice(umis[:6]), random.choice(umis[:6])
        else:  # Remaining reads use all UMIs
            umi1, umi2 = random.choice(umis), random.choice(umis)
        
        umi_distribution[umi1] += 1
        umi_distribution[umi2] += 1
        
        # Build read by replacing N's with actual UMIs and adding mutations
        read_seq = (
            add_mutation(prefix_seq) + 
            umi1 + 
            add_mutation(middle_seq) + 
            umi2 + 
            add_mutation(suffix_seq)
        )
        
        read_qual = generate_quality_scores(len(read_seq))
        reads.append(f"@read_{read_id:04d}\n{read_seq}\n+\n{read_qual}")
    
    # Write FASTQ
    fastq_path = Path("data/1-fastq_files/demo_sample.fastq")
    fastq_path.write_text("\n".join(reads) + "\n")
    print(f"✅ FASTQ created: {fastq_path} ({len(reads)} reads)")
    
    # Print UMI distribution
    print(f"\n📊 UMI distribution (top 5):")
    for umi, count in sorted(umi_distribution.items(), key=lambda x: -x[1])[:5]:
        print(f"   {umi}: {count} occurrences")
    
    print(f"\n✅ All synthetic demo data generated successfully!")
    print(f"   Total size: ~{len(reads) * len(reads[0]) / 1024:.1f} KB")
    print(f"\n   Reference structure:")
    print(f"   - Prefix: 0-{prefix_len-1}")
    print(f"   - UMI1: {prefix_len}-{prefix_len+11}")
    print(f"   - Middle: {prefix_len+12}-{prefix_len+12+middle_len-1}")
    print(f"   - UMI2: {prefix_len+12+middle_len}-{prefix_len+12+middle_len+11}")
    print(f"   - Suffix: {prefix_len+12+middle_len+12}-{len(ref_seq)-1}")

if __name__ == "__main__":
    main()

