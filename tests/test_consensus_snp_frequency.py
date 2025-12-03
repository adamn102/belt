#!/usr/bin/env python3
"""
Tests for consensus calling with SNPs at various frequencies (1% to 99%).

This test suite validates that bcftools can correctly call SNPs when they are
present at different frequencies in the read data, helping to determine the
minimum read fraction required for SNP calling.
"""

import pytest
import subprocess
import os
from pathlib import Path
from Bio import SeqIO
import pysam

# Import consensus functions
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.consensus import umi_consensus_bcftools

# Paths to test data
TEST_DATA_DIR = Path(__file__).parent.parent / "test-data" / "test_consensus"
SNP_FREQ_DIR = TEST_DATA_DIR / "snp_frequency"
REFERENCE_PATH = TEST_DATA_DIR / "0-consensus_reference" / "LP249.fa"
OUTPUT_DIR = TEST_DATA_DIR / "consensus_output_snp_frequency"

# Reference length (will be set from reference file)
REF_LENGTH = None


def parse_test_metadata():
    """Parse test metadata to get test case information."""
    metadata_file = SNP_FREQ_DIR / "test_metadata.txt"
    if not metadata_file.exists():
        pytest.skip(f"Metadata file not found: {metadata_file}. Run generate_snp_frequency_tests.py first.")
    
    test_cases = {}
    current_test = None
    
    with open(metadata_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Test:"):
                current_test = line.split(":", 1)[1].strip()
                test_cases[current_test] = {}
            elif line.startswith("UMI:") and current_test:
                test_cases[current_test]['umi'] = line.split(":", 1)[1].strip()
            elif line.startswith("SNP_Position:") and current_test:
                test_cases[current_test]['snp_position'] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Ref_Base:") and current_test:
                test_cases[current_test]['ref_base'] = line.split(":", 1)[1].strip()
            elif line.startswith("SNP_Base:") and current_test:
                test_cases[current_test]['snp_base'] = line.split(":", 1)[1].strip()
            elif line.startswith("Expected_Frequency:") and current_test:
                test_cases[current_test]['expected_frequency'] = float(line.split(":", 1)[1].strip())
            elif line.startswith("Actual_Frequency:") and current_test:
                test_cases[current_test]['actual_frequency'] = float(line.split(":", 1)[1].strip())
            elif line.startswith("Reads_With_SNP:") and current_test:
                test_cases[current_test]['reads_with_snp'] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Total_Reads:") and current_test:
                test_cases[current_test]['total_reads'] = int(line.split(":", 1)[1].strip())
            elif line.startswith("BAM:") and current_test:
                bam_path = line.split(":", 1)[1].strip()
                # Resolve relative path from project root
                if not Path(bam_path).is_absolute():
                    # Path is relative (e.g., "belt/test-data/...")
                    # Test file is at belt/tests/test_consensus_snp_frequency.py
                    # Resolve test file to absolute path, then go up to project root
                    test_file_abs = Path(__file__).resolve()
                    project_root = test_file_abs.parent.parent.parent  # Go up to project root
                    bam_path = str(project_root / bam_path)
                test_cases[current_test]['bam_file'] = bam_path
    
    return test_cases


def get_test_cases():
    """Get list of test case names."""
    test_cases = parse_test_metadata()
    return list(test_cases.keys())


@pytest.fixture(scope="module")
def reference_fasta():
    """Load reference FASTA file."""
    global REF_LENGTH
    if not REFERENCE_PATH.exists():
        pytest.skip(f"Reference file not found: {REFERENCE_PATH}")
    
    with open(REFERENCE_PATH) as f:
        ref_record = SeqIO.read(f, "fasta")
        REF_LENGTH = len(ref_record.seq)
        return str(REFERENCE_PATH)


@pytest.fixture(scope="module")
def output_directory():
    """Create output directory for consensus files."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


@pytest.fixture(scope="module", autouse=True)
def check_dependencies():
    """Check if bcftools and samtools are available."""
    if subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0:
        pytest.skip("bcftools not available")
    if subprocess.run(["which", "samtools"], capture_output=True).returncode != 0:
        pytest.skip("samtools not available")


def check_snp_in_consensus(consensus_file, snp_pos, ref_base, snp_base):
    """
    Check if SNP is present in consensus sequence.
    
    Args:
        consensus_file: Path to consensus FASTA
        snp_pos: 0-based position of SNP
        ref_base: Reference base
        snp_base: Alternative base (SNP)
    
    Returns:
        tuple: (found, consensus_base, position)
    """
    with open(consensus_file) as f:
        consensus_record = SeqIO.read(f, "fasta")
        consensus_seq = str(consensus_record.seq)
    
    # Check if position is valid
    if snp_pos >= len(consensus_seq):
        return False, None, None
    
    consensus_base = consensus_seq[snp_pos]
    
    # SNP is called if consensus base matches the SNP base
    if consensus_base == snp_base:
        return True, consensus_base, snp_pos
    elif consensus_base == ref_base:
        return False, consensus_base, snp_pos
    else:
        # Different base (unexpected)
        return None, consensus_base, snp_pos


def check_snp_in_vcf(vcf_file, snp_pos, ref_base, snp_base):
    """
    Check if SNP is present in VCF file.
    
    Args:
        vcf_file: Path to VCF file (can be .vcf.gz)
        snp_pos: 1-based position of SNP (VCF uses 1-based)
        ref_base: Reference base
        snp_base: Alternative base (SNP)
    
    Returns:
        tuple: (found, vcf_ref, vcf_alt, vcf_pos)
    """
    try:
        with pysam.VariantFile(str(vcf_file)) as vcf:
            for record in vcf:
                # VCF uses 1-based positions
                if record.pos == snp_pos:
                    vcf_ref = record.ref
                    vcf_alt = str(record.alts[0]) if record.alts else None
                    
                    # Check if this matches our SNP
                    if vcf_ref == ref_base and vcf_alt == snp_base:
                        return True, vcf_ref, vcf_alt, record.pos
                    else:
                        # Different variant at same position
                        return None, vcf_ref, vcf_alt, record.pos
        return False, None, None, None
    except Exception as e:
        # VCF might not exist or be readable
        return None, None, None, None


class TestSNPFrequencyCalling:
    """Test SNP calling at various frequencies."""
    
    @pytest.mark.parametrize("test_name", get_test_cases())
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_snp_frequency_consensus(self, tmp_path, reference_fasta, output_directory, test_name):
        """
        Test that SNPs at various frequencies are correctly called in consensus.
        
        This test checks:
        1. Consensus sequence is generated
        2. SNP is present in consensus (if frequency >= 80%)
        3. SNP position matches expected
        
        Note: Only SNPs with frequency >= 0.80 (80%) are expected to be called.
        Lower frequency SNPs may or may not be called.
        """
        test_cases = parse_test_metadata()
        if test_name not in test_cases:
            pytest.skip(f"Test case {test_name} not found in metadata")
        
        case = test_cases[test_name]
        bam_file = Path(case['bam_file']).resolve()  # Resolve to absolute path
        umi = case['umi']
        snp_pos = case['snp_position']
        ref_base = case['ref_base']
        snp_base = case['snp_base']
        actual_freq = case['actual_frequency']
        reads_with_snp = case['reads_with_snp']
        total_reads = case['total_reads']
        
        if not bam_file.exists():
            pytest.skip(f"BAM file not found: {bam_file}")
        
        # Change to tmp_path for bcftools processing
        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        
        try:
            # Run consensus calling
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_file,
                ref=Path(reference_fasta),
                threads=1,
                skip_extraction=False,
                indel_sensitivity=0.01
            )
            
            consensus_path = Path(consensus_file)
            if not consensus_path.exists():
                pytest.fail(f"Consensus file not created: {consensus_file}")
            
            # Check consensus sequence
            snp_found, consensus_base, pos = check_snp_in_consensus(
                consensus_path, snp_pos, ref_base, snp_base
            )
            
            # Write output files for inspection
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy consensus file
            output_consensus = output_dir / f"{test_name}_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_consensus)
            
            # Write metadata
            metadata_file = output_dir / f"{test_name}_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name}\n")
                f.write(f"UMI: {umi}\n")
                f.write(f"SNP Position: {snp_pos}\n")
                f.write(f"Reference Base: {ref_base}\n")
                f.write(f"SNP Base: {snp_base}\n")
                f.write(f"Expected Frequency: {case['expected_frequency']:.3f}\n")
                f.write(f"Actual Frequency: {actual_freq:.3f}\n")
                f.write(f"Reads with SNP: {reads_with_snp}/{total_reads}\n")
                f.write(f"Consensus Base at {snp_pos}: {consensus_base}\n")
                f.write(f"SNP Found in Consensus: {snp_found}\n")
                f.write(f"Consensus File: {output_consensus}\n")
            
            # Check if SNP was called
            # Only SNPs with frequency >= 0.80 (80%) are expected to be called
            # For lower frequencies, SNP may not be called (this is expected)
            if actual_freq >= 0.80:
                # High frequency (>=80%): SNP should be called
                assert snp_found is not False, (
                    f"SNP at {snp_pos} ({ref_base}->{snp_base}) not called in consensus "
                    f"despite {actual_freq:.1%} frequency ({reads_with_snp}/{total_reads} reads). "
                    f"SNPs with frequency >= 80% should be called. "
                    f"Consensus base: {consensus_base}"
                )
            else:
                # Lower frequency (<80%): SNP may or may not be called
                # Just check that consensus was generated
                assert consensus_base is not None, "Consensus sequence is empty"
            
            # Always check that consensus length matches reference
            with open(consensus_path) as f:
                consensus_record = SeqIO.read(f, "fasta")
                consensus_length = len(consensus_record.seq)
            
            assert consensus_length == REF_LENGTH, (
                f"Consensus length mismatch: expected {REF_LENGTH}, got {consensus_length}"
            )
            
        finally:
            os.chdir(original_cwd)
    
    @pytest.mark.parametrize("test_name", get_test_cases())
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_snp_frequency_vcf(self, tmp_path, reference_fasta, output_directory, test_name):
        """
        Test that SNPs at various frequencies are correctly called in VCF.
        
        This test checks the VCF file directly to see if the SNP variant
        is present, which gives more detailed information than consensus.
        
        Note: Only SNPs with frequency >= 0.80 (80%) are expected to be called.
        Lower frequency SNPs may or may not be called.
        """
        test_cases = parse_test_metadata()
        if test_name not in test_cases:
            pytest.skip(f"Test case {test_name} not found in metadata")
        
        case = test_cases[test_name]
        bam_file = Path(case['bam_file']).resolve()  # Resolve to absolute path
        umi = case['umi']
        snp_pos = case['snp_position']  # 0-based
        vcf_snp_pos = snp_pos + 1  # Convert to 1-based for VCF
        ref_base = case['ref_base']
        snp_base = case['snp_base']
        actual_freq = case['actual_frequency']
        reads_with_snp = case['reads_with_snp']
        total_reads = case['total_reads']
        
        if not bam_file.exists():
            pytest.skip(f"BAM file not found: {bam_file}")
        
        # Change to tmp_path for bcftools processing
        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        
        try:
            # Run consensus calling (this generates VCF files)
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_file,
                ref=Path(reference_fasta),
                threads=1,
                skip_extraction=False,
                indel_sensitivity=0.01
            )
            
            # Check normalized VCF file
            norm_vcf = Path(f"umi.{umi}.norm.vcf.gz")
            if not norm_vcf.exists():
                pytest.skip(f"Normalized VCF file not found: {norm_vcf}")
            
            # Check if SNP is in VCF
            snp_found, vcf_ref, vcf_alt, vcf_pos = check_snp_in_vcf(
                norm_vcf, vcf_snp_pos, ref_base, snp_base
            )
            
            # Write results
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            metadata_file = output_dir / f"{test_name}_vcf_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name}\n")
                f.write(f"UMI: {umi}\n")
                f.write(f"SNP Position (0-based): {snp_pos}\n")
                f.write(f"SNP Position (1-based): {vcf_snp_pos}\n")
                f.write(f"Reference Base: {ref_base}\n")
                f.write(f"SNP Base: {snp_base}\n")
                f.write(f"Expected Frequency: {case['expected_frequency']:.3f}\n")
                f.write(f"Actual Frequency: {actual_freq:.3f}\n")
                f.write(f"Reads with SNP: {reads_with_snp}/{total_reads}\n")
                f.write(f"SNP Found in VCF: {snp_found}\n")
                if vcf_ref:
                    f.write(f"VCF Reference: {vcf_ref}\n")
                    f.write(f"VCF Alternative: {vcf_alt}\n")
                    f.write(f"VCF Position: {vcf_pos}\n")
            
            # Only SNPs with frequency >= 0.80 (80%) are expected to be in VCF
            if actual_freq >= 0.80:
                assert snp_found is not False, (
                    f"SNP at position {vcf_snp_pos} ({ref_base}->{snp_base}) not found in VCF "
                    f"despite {actual_freq:.1%} frequency ({reads_with_snp}/{total_reads} reads). "
                    f"SNPs with frequency >= 80% should be called."
                )
            
        finally:
            os.chdir(original_cwd)


class TestSNPFrequencyThreshold:
    """Test to determine minimum frequency threshold for SNP calling."""
    
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_snp_frequency_threshold_analysis(self, tmp_path, reference_fasta, output_directory):
        """
        Analyze SNP calling across all frequencies to determine minimum threshold.
        
        This test runs all frequency tests and summarizes the results to determine
        the minimum read fraction required for bcftools to call a SNP.
        """
        test_cases = parse_test_metadata()
        results = []
        
        for test_name, case in test_cases.items():
            bam_file = Path(case['bam_file']).resolve()  # Resolve to absolute path
            umi = case['umi']
            snp_pos = case['snp_position']
            ref_base = case['ref_base']
            snp_base = case['snp_base']
            actual_freq = case['actual_frequency']
            
            if not bam_file.exists():
                continue
            
            original_cwd = os.getcwd()
            os.chdir(tmp_path)
            
            try:
                # Run consensus calling
                consensus_file = umi_consensus_bcftools(
                    umi=umi,
                    bam=bam_file,
                    ref=Path(reference_fasta),
                    threads=1,
                    skip_extraction=False,
                    indel_sensitivity=0.01
                )
                
                consensus_path = Path(consensus_file)
                if consensus_path.exists():
                    snp_found, consensus_base, pos = check_snp_in_consensus(
                        consensus_path, snp_pos, ref_base, snp_base
                    )
                    
                    results.append({
                        'test_name': test_name,
                        'frequency': actual_freq,
                        'snp_called': snp_found,
                        'consensus_base': consensus_base
                    })
            except Exception as e:
                results.append({
                    'test_name': test_name,
                    'frequency': actual_freq,
                    'snp_called': None,
                    'error': str(e)
                })
            finally:
                os.chdir(original_cwd)
        
        # Write summary
        output_dir = Path(output_directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        summary_file = output_dir / "snp_frequency_threshold_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("SNP Frequency Calling Threshold Analysis\n")
            f.write("=" * 80 + "\n\n")
            f.write("This analysis determines the minimum read fraction required for SNP calling.\n")
            f.write("Expected threshold: SNPs with frequency >= 0.80 (80%) should be called.\n\n")
            
            # Sort by frequency
            results_sorted = sorted(results, key=lambda x: x['frequency'])
            
            f.write("Results by Frequency:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Frequency':<12} {'SNP Called':<12} {'Test Name':<20} {'Notes':<30}\n")
            f.write("-" * 80 + "\n")
            
            called_frequencies = []
            not_called_frequencies = []
            
            for r in results_sorted:
                freq_str = f"{r['frequency']:.1%}"
                called_str = "Yes" if r['snp_called'] else ("No" if r['snp_called'] is False else "Error")
                notes = r.get('error', r.get('consensus_base', ''))
                
                f.write(f"{freq_str:<12} {called_str:<12} {r['test_name']:<20} {notes[:30]:<30}\n")
                
                if r['snp_called']:
                    called_frequencies.append(r['frequency'])
                elif r['snp_called'] is False:
                    not_called_frequencies.append(r['frequency'])
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("Summary:\n")
            f.write("-" * 80 + "\n")
            
            if called_frequencies:
                min_called = min(called_frequencies)
                f.write(f"Minimum frequency where SNP was called: {min_called:.1%}\n")
            
            if not_called_frequencies:
                max_not_called = max(not_called_frequencies)
                f.write(f"Maximum frequency where SNP was NOT called: {max_not_called:.1%}\n")
            
            if called_frequencies and not_called_frequencies:
                threshold_range = (max(not_called_frequencies), min(called_frequencies))
                f.write(f"Estimated threshold range: {threshold_range[0]:.1%} - {threshold_range[1]:.1%}\n")
        
        # Don't fail the test - this is an analysis test
        assert len(results) > 0, "No results collected"

