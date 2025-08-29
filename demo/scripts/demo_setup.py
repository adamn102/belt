#!/usr/bin/env python3
"""
BELT Demo Setup Script
This script prepares the demo environment by setting up directories and copying test data.
"""

import os
import shutil
import sys
from pathlib import Path

def main():
    """Setup the demo environment"""
    
    print("🧪 Setting up BELT Demo Pipeline")
    print("=" * 50)
    
    # Get the demo directory
    demo_dir = Path(__file__).parent.parent
    project_root = demo_dir.parent
    
    # Create data directories
    data_dirs = [
        "data/0-reference_sequence",
        "data/0-library_whitelist", 
        "data/1-fastq_files",
        "data/2-fastq_filtered",
        "data/3-mapped_reads",
        "data/4-sorted_reads",
        "data/5-umi_maps",
        "data/6-umi_groups",
        "data/7-consensus_calls",
        "data/8-compiled_consensus",
        "logs/fastq",
        "logs/minimap2",
        "logs/consensus"
    ]
    
    print("📁 Creating data directories...")
    for dir_path in data_dirs:
        full_path = demo_dir / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"   ✅ {dir_path}")
    
    # Copy test data
    print("\n📋 Copying test data...")
    
    # Copy reference sequence
    ref_source = project_root / "scripts" / "test-data" / "test-umi_map" / "reference_sequence" / "B145.fa"
    ref_dest = demo_dir / "data" / "0-reference_sequence" / "B145-XTMMF9-12.fa"
    if ref_source.exists():
        shutil.copy2(ref_source, ref_dest)
        print(f"   ✅ Reference: {ref_source.name}")
    else:
        print(f"   ❌ Reference not found: {ref_source}")
        return
    
    # Copy BAM file (we'll need to create a FASTQ from it)
    bam_source = project_root / "scripts" / "test-data" / "test-umi_map" / "sorted_reads" / "B145-XTMMF9-12.bam"
    fastq_dest = demo_dir / "data" / "1-fastq_files" / "B145-XTMMF9-12.fastq"
    
    if bam_source.exists():
        print(f"   ✅ BAM file: {bam_source.name}")
        print("   🔄 Converting BAM to FASTQ...")
        
        # Convert BAM to FASTQ using samtools
        try:
            import subprocess
            cmd = f"samtools fastq {bam_source} > {fastq_dest}"
            subprocess.run(cmd, shell=True, check=True)
            print(f"   ✅ FASTQ created: {fastq_dest.name}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("   ⚠️  samtools not available, creating dummy FASTQ...")
            # Create a dummy FASTQ file for demo purposes
            create_dummy_fastq(fastq_dest)
    else:
        print(f"   ❌ BAM file not found: {bam_source}")
        print("   🔄 Creating dummy FASTQ...")
        create_dummy_fastq(fastq_dest)
    
    print("\n✅ Demo setup complete!")
    print("\n📝 Next steps:")
    print("   1. Install required tools: minimap2, samtools, fastplong")
    print("   2. Run the pipeline: snakemake --cores 4")
    print("   3. Check results in data/ directories")

def create_dummy_fastq(fastq_path):
    """Create a dummy FASTQ file for demo purposes"""
    dummy_reads = [
        ("@read_1", "ATCGATCGATCG", "+", "IIIIIIIIIIII"),
        ("@read_2", "GCTAGCTAGCTA", "+", "IIIIIIIIIIII"),
        ("@read_3", "TAGCTAGCTAGC", "+", "IIIIIIIIIIII")
    ]
    
    with open(fastq_path, 'w') as f:
        for read_id, seq, plus, qual in dummy_reads:
            f.write(f"{read_id}\n{seq}\n{plus}\n{qual}\n")
    
    print(f"   ✅ Dummy FASTQ created: {fastq_path.name}")

if __name__ == "__main__":
    main()
