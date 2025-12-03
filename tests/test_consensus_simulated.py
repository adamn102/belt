#!/usr/bin/env python3
"""
Tests for consensus calling using simulated BAM files.

This test suite validates consensus calling on comprehensive simulated data
with known mutations (SNPs, insertions, deletions, combinations).
"""

import pytest
import subprocess
import tempfile
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import consensus functions
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.consensus import umi_consensus_bcftools
from scripts.consensus_improved import umi_consensus_with_prefiltering, hybrid_consensus_calling


# Paths to test data
TEST_DATA_DIR = Path(__file__).parent.parent / "test-data" / "test_consensus"
SIMULATED_DIR = TEST_DATA_DIR / "simulated"
REFERENCE_PATH = TEST_DATA_DIR / "0-consensus_reference" / "LP249.fa"
OUTPUT_DIR = TEST_DATA_DIR / "consensus_output_simulated"


def parse_test_metadata():
    """Parse test metadata to get test case information."""
    metadata_file = SIMULATED_DIR / "test_metadata.txt"
    if not metadata_file.exists():
        pytest.skip(f"Metadata file not found: {metadata_file}")
    
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
            elif line.startswith("Description:") and current_test:
                test_cases[current_test]['description'] = line.split(":", 1)[1].strip()
            elif line.startswith("Expected:") and current_test:
                test_cases[current_test]['expected'] = line.split(":", 1)[1].strip()
            elif line.startswith("Reads:") and current_test:
                test_cases[current_test]['reads'] = int(line.split(":", 1)[1].strip())
            elif line.startswith("BAM:") and current_test:
                test_cases[current_test]['bam'] = Path(__file__).parent.parent.parent / line.split(":", 1)[1].strip()
    
    return test_cases


def get_reference_sequence():
    """Get reference sequence and length."""
    with open(REFERENCE_PATH) as f:
        ref_record = SeqIO.read(f, "fasta")
        return str(ref_record.seq), len(ref_record.seq)


# Parse test cases once
TEST_CASES = parse_test_metadata()
REF_SEQ, REF_LENGTH = get_reference_sequence()


@pytest.fixture(scope="module")
def output_directory():
    """Create output directory for consensus files."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


@pytest.fixture(scope="module")
def reference_fasta():
    """Reference FASTA file."""
    assert REFERENCE_PATH.exists(), f"Reference not found: {REFERENCE_PATH}"
    return REFERENCE_PATH


class TestSimulatedSNPs:
    """Test consensus calling on simulated SNP data."""
    
    @pytest.mark.parametrize("test_name", [
        "SNP_1",
        "SNP_5", 
        "SNP_20"
    ])
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_snp_consensus(self, tmp_path, reference_fasta, output_directory, test_name):
        """Test consensus calling on SNP test cases."""
        if test_name not in TEST_CASES:
            pytest.skip(f"Test case {test_name} not found")
        
        test_case = TEST_CASES[test_name]
        bam_path = test_case['bam']
        umi = test_case['umi']
        
        assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            
            # Run consensus calling
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_path.resolve(),
                ref=reference_fasta.resolve(),
                threads=1,
                indel_sensitivity=0.1
            )
            
            consensus_path = Path(consensus_file)
            assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
            # Read consensus
            with open(consensus_path) as f:
                consensus_lines = f.readlines()
                consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
            consensus_length = len(consensus_seq)
            
            # For SNPs, consensus should be same length as reference
            assert consensus_length == REF_LENGTH, \
                f"SNP consensus length mismatch: expected {REF_LENGTH}, got {consensus_length}"
            
            # Write to output directory for inspection
            output_file = output_directory / f"{test_name}_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_file)
            
            # Write metadata
            metadata_file = output_directory / f"{test_name}_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name}\n")
                f.write(f"Description: {test_case['description']}\n")
                f.write(f"Expected: {test_case['expected']}\n")
                f.write(f"Reference length: {REF_LENGTH} bp\n")
                f.write(f"Consensus length: {consensus_length} bp\n")
                f.write(f"Length match: {'✓' if consensus_length == REF_LENGTH else '✗'}\n")
                f.write(f"\nConsensus file: {output_file.name}\n")
            
            print(f"\n✓ {test_name}: Consensus length {consensus_length} bp (expected {REF_LENGTH} bp)")
            
        finally:
            os.chdir(original_cwd)


class TestSimulatedInsertions:
    """Test consensus calling on simulated insertion data."""
    
    @pytest.mark.parametrize("test_name", [
        "INS_1", "INS_2", "INS_5", "INS_10", 
        "INS_25", "INS_50", "INS_100", "INS_200"
    ])
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_insertion_consensus(self, tmp_path, reference_fasta, output_directory, test_name):
        """Test consensus calling on insertion test cases."""
        if test_name not in TEST_CASES:
            pytest.skip(f"Test case {test_name} not found")
        
        test_case = TEST_CASES[test_name]
        bam_path = test_case['bam']
        umi = test_case['umi']
        
        # Extract expected insertion size from test name
        expected_ins_size = int(test_name.split('_')[1])
        
        assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            
            # Run consensus calling
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_path.resolve(),
                ref=reference_fasta.resolve(),
                threads=1,
                indel_sensitivity=0.1
            )
            
            consensus_path = Path(consensus_file)
            assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
            # Read consensus
            with open(consensus_path) as f:
                consensus_lines = f.readlines()
                consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
            consensus_length = len(consensus_seq)
            length_diff = consensus_length - REF_LENGTH
            
            # For insertions, consensus should be longer than reference
            assert consensus_length >= REF_LENGTH, \
                f"Insertion consensus should not be shorter: {consensus_length} < {REF_LENGTH}"
            
            # Strict validation: length difference must match expected (within 2bp tolerance)
            expected_length_diff = expected_ins_size
            length_match = abs(length_diff - expected_length_diff) <= 2
            
            # Write to output directory
            output_file = output_directory / f"{test_name}_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_file)
            
            # Write metadata
            metadata_file = output_directory / f"{test_name}_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name}\n")
                f.write(f"Description: {test_case['description']}\n")
                f.write(f"Expected insertion size: {expected_ins_size} bp\n")
                f.write(f"Reference length: {REF_LENGTH} bp\n")
                f.write(f"Consensus length: {consensus_length} bp\n")
                f.write(f"Length difference: {length_diff} bp\n")
                f.write(f"Expected difference: {expected_length_diff} bp\n")
                f.write(f"Match: {'✓ PASS' if length_match else '✗ FAIL'}\n")
                f.write(f"\nConsensus file: {output_file.name}\n")
            
            print(f"\n{test_name}: Consensus {consensus_length} bp (+{length_diff} bp, expected +{expected_length_diff} bp) - {'✓ PASS' if length_match else '✗ FAIL'}")
            
            # Assert strict match
            assert length_match, \
                f"Insertion not correctly called: expected length difference {expected_length_diff} bp, got {length_diff} bp"
            
        finally:
            os.chdir(original_cwd)


class TestSimulatedDeletions:
    """Test consensus calling on simulated deletion data."""
    
    @pytest.mark.parametrize("test_name", [
        "DEL_1", "DEL_2", "DEL_5", "DEL_10", "DEL_25", "DEL_50",
        "DEL_100", "DEL_200", "DEL_300", "DEL_400", "DEL_500"
    ])
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_deletion_consensus(self, tmp_path, reference_fasta, output_directory, test_name):
        """Test consensus calling on deletion test cases."""
        if test_name not in TEST_CASES:
            pytest.skip(f"Test case {test_name} not found")
        
        test_case = TEST_CASES[test_name]
        bam_path = test_case['bam']
        umi = test_case['umi']
        
        # Extract expected deletion size from test name
        expected_del_size = int(test_name.split('_')[1])
        
        assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            
            # Run consensus calling
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_path.resolve(),
                ref=reference_fasta.resolve(),
                threads=1,
                indel_sensitivity=0.1
            )
            
            consensus_path = Path(consensus_file)
            assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
            # Read consensus
            with open(consensus_path) as f:
                consensus_lines = f.readlines()
                consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
            consensus_length = len(consensus_seq)
            length_diff = REF_LENGTH - consensus_length
            
            # For deletions, consensus should be shorter than reference
            assert consensus_length <= REF_LENGTH, \
                f"Deletion consensus should not be longer: {consensus_length} > {REF_LENGTH}"
            
            # Strict validation: length difference must match expected (within 2bp tolerance)
            expected_length_diff = expected_del_size
            length_match = abs(length_diff - expected_length_diff) <= 2
            
            # Write to output directory
            output_file = output_directory / f"{test_name}_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_file)
            
            # Write metadata
            metadata_file = output_directory / f"{test_name}_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name}\n")
                f.write(f"Description: {test_case['description']}\n")
                f.write(f"Expected deletion size: {expected_del_size} bp\n")
                f.write(f"Reference length: {REF_LENGTH} bp\n")
                f.write(f"Consensus length: {consensus_length} bp\n")
                f.write(f"Length difference: {length_diff} bp\n")
                f.write(f"Expected difference: {expected_length_diff} bp\n")
                f.write(f"Match: {'✓ PASS' if length_match else '✗ FAIL'}\n")
                f.write(f"\nConsensus file: {output_file.name}\n")
            
            print(f"\n{test_name}: Consensus {consensus_length} bp (-{length_diff} bp, expected -{expected_length_diff} bp) - {'✓ PASS' if length_match else '✗ FAIL'}")
            
            # Assert strict match
            assert length_match, \
                f"Deletion not correctly called: expected length difference {expected_length_diff} bp, got {length_diff} bp"
            
        finally:
            os.chdir(original_cwd)


class TestSizeVariableDeletion:
    """Test consensus calling on size-variable deletion (challenging case)."""
    
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_size_variable_deletion_standard(self, tmp_path, reference_fasta, output_directory):
        """Test standard consensus calling on size-variable deletion."""
        test_name = "DEL_VARIABLE"
        if test_name not in TEST_CASES:
            pytest.skip(f"Test case {test_name} not found")
        
        test_case = TEST_CASES[test_name]
        bam_path = test_case['bam']
        umi = test_case['umi']
        
        assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            
            # Run standard consensus calling
            consensus_file = umi_consensus_bcftools(
                umi=umi,
                bam=bam_path.resolve(),
                ref=reference_fasta.resolve(),
                threads=1,
                indel_sensitivity=0.05  # Lower for better detection
            )
            
            consensus_path = Path(consensus_file)
            assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
            # Read consensus
            with open(consensus_path) as f:
                consensus_lines = f.readlines()
                consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
            consensus_length = len(consensus_seq)
            length_diff = REF_LENGTH - consensus_length
            
            # Write to output directory
            output_file = output_directory / f"{test_name}_standard_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_file)
            
            # Write metadata
            metadata_file = output_directory / f"{test_name}_standard_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name} (Standard Method)\n")
                f.write(f"Description: {test_case['description']}\n")
                f.write(f"Expected: {test_case['expected']}\n")
                f.write(f"Reference length: {REF_LENGTH} bp\n")
                f.write(f"Consensus length: {consensus_length} bp\n")
                f.write(f"Length difference: {length_diff} bp\n")
                f.write(f"Expected: ~400-500bp deletion (most common: 500bp)\n")
                if length_diff > 300:
                    f.write(f"Result: ✓ Large deletion appears to be called\n")
                elif length_diff > 50:
                    f.write(f"Result: ? Partial deletion may be called\n")
                else:
                    f.write(f"Result: ⚠ Deletion may not be called (size variation challenge)\n")
                f.write(f"\nConsensus file: {output_file.name}\n")
            
            print(f"\n{test_name} (standard): Consensus {consensus_length} bp (-{length_diff} bp)")
            
        finally:
            os.chdir(original_cwd)
    
    @pytest.mark.skipif(
        subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
        reason="bcftools not available"
    )
    def test_size_variable_deletion_prefiltered(self, tmp_path, reference_fasta, output_directory):
        """Test pre-filtered consensus calling on size-variable deletion."""
        test_name = "DEL_VARIABLE"
        if test_name not in TEST_CASES:
            pytest.skip(f"Test case {test_name} not found")
        
        test_case = TEST_CASES[test_name]
        bam_path = test_case['bam']
        umi = test_case['umi']
        
        assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
        original_cwd = os.getcwd()
        try:
            os.chdir(tmp_path)
            
            # Run pre-filtered consensus calling
            consensus_file = umi_consensus_with_prefiltering(
                umi=umi,
                bam=bam_path.resolve(),
                ref=reference_fasta.resolve(),
                threads=1,
                deletion_size_tolerance=50,
                use_filtered=True,
                indel_sensitivity=0.05
            )
            
            consensus_path = Path(consensus_file)
            assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
            # Read consensus
            with open(consensus_path) as f:
                consensus_lines = f.readlines()
                consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
            consensus_length = len(consensus_seq)
            length_diff = REF_LENGTH - consensus_length
            
            # Write to output directory
            output_file = output_directory / f"{test_name}_prefiltered_consensus.fasta"
            import shutil
            shutil.copy2(consensus_path, output_file)
            
            # Write metadata
            metadata_file = output_directory / f"{test_name}_prefiltered_metadata.txt"
            with open(metadata_file, 'w') as f:
                f.write(f"Test: {test_name} (Pre-filtered Method)\n")
                f.write(f"Description: {test_case['description']}\n")
                f.write(f"Expected: {test_case['expected']}\n")
                f.write(f"Reference length: {REF_LENGTH} bp\n")
                f.write(f"Consensus length: {consensus_length} bp\n")
                f.write(f"Length difference: {length_diff} bp\n")
                f.write(f"Expected: ~400-500bp deletion (most common: 500bp)\n")
                if length_diff > 300:
                    f.write(f"Result: ✓ Large deletion appears to be called\n")
                elif length_diff > 50:
                    f.write(f"Result: ? Partial deletion may be called\n")
                else:
                    f.write(f"Result: ⚠ Deletion may not be called\n")
                f.write(f"\nConsensus file: {output_file.name}\n")
            
            print(f"\n{test_name} (pre-filtered): Consensus {consensus_length} bp (-{length_diff} bp)")
            
        finally:
            os.chdir(original_cwd)


# class TestCombinedMutations:
#     """Test consensus calling on combined mutation test cases."""
    
#     def get_expected_mutations(self, test_name):
#         """Parse expected mutations from test case description."""
#         if test_name == "COMBINED_SNP_DEL":
#             # 3 SNPs + 50bp deletion
#             # Expected: length_diff = -50bp (deletion)
#             # Expected: SNPs at positions 300, 600, 900 (from generation script)
#             return {
#                 'expected_length_diff': -50,
#                 'expected_snp_positions': [300, 600, 900],  # Approximate, from generation script
#                 'expected_deletion_size': 50
#             }
#         elif test_name == "COMBINED_INDELS":
#             # 8bp insertion + 20bp deletion
#             # Expected: length_diff = +8 - 20 = -12bp (net deletion)
#             return {
#                 'expected_length_diff': -12,  # 8bp insertion - 20bp deletion
#                 'expected_insertion_size': 8,
#                 'expected_deletion_size': 20
#             }
#         return {}
    
#     def find_snps_in_consensus(self, ref_seq, consensus_seq):
#         """Find SNP positions by comparing reference and consensus sequences."""
#         snps = []
#         min_len = min(len(ref_seq), len(consensus_seq))
#         for i in range(min_len):
#             if ref_seq[i] != consensus_seq[i]:
#                 snps.append(i)
#         return snps
    
#     @pytest.mark.parametrize("test_name", [
#         "COMBINED_SNP_DEL",
#         "COMBINED_INDELS"
#     ])
#     @pytest.mark.skipif(
#         subprocess.run(["which", "bcftools"], capture_output=True).returncode != 0,
#         reason="bcftools not available"
#     )
#     def test_combined_mutations(self, tmp_path, reference_fasta, output_directory, test_name):
#         """Test consensus calling on combined mutation test cases.
        
#         Tests should only PASS if:
#         - Expected length difference matches actual length difference AND/OR
#         - Expected SNP positions match actual SNP positions (for SNP-containing tests)
#         """
#         if test_name not in TEST_CASES:
#             pytest.skip(f"Test case {test_name} not found")
        
#         test_case = TEST_CASES[test_name]
#         bam_path = test_case['bam']
#         umi = test_case['umi']
#         expected = self.get_expected_mutations(test_name)
        
#         assert bam_path.exists(), f"BAM file not found: {bam_path}"
        
#         # Read reference sequence
#         with open(reference_fasta) as f:
#             ref_record = SeqIO.read(f, "fasta")
#             ref_seq = str(ref_record.seq)
        
#         original_cwd = os.getcwd()
#         try:
#             os.chdir(tmp_path)
            
#             # Run consensus calling
#             consensus_file = umi_consensus_bcftools(
#                 umi=umi,
#                 bam=bam_path.resolve(),
#                 ref=reference_fasta.resolve(),
#                 threads=1,
#                 indel_sensitivity=0.1
#             )
            
#             consensus_path = Path(consensus_file)
#             assert consensus_path.exists(), f"Consensus file not created for {test_name}"
            
#             # Read consensus
#             with open(consensus_path) as f:
#                 consensus_lines = f.readlines()
#                 consensus_seq = ''.join(line.strip() for line in consensus_lines[1:])
            
#             consensus_length = len(consensus_seq)
#             length_diff = consensus_length - REF_LENGTH
            
#             # Write to output directory
#             output_file = output_directory / f"{test_name}_consensus.fasta"
#             import shutil
#             shutil.copy2(consensus_path, output_file)
            
#             # Validate results based on test type
#             passed = False
#             validation_details = []
            
#             if test_name == "COMBINED_SNP_DEL":
#                 # Check length difference (should be -50bp due to deletion)
#                 length_match = abs(length_diff - expected['expected_length_diff']) <= 2  # Allow 2bp tolerance
#                 validation_details.append(f"Length difference: {length_diff} bp (expected: {expected['expected_length_diff']} bp) - {'✓' if length_match else '✗'}")
                
#                 # Check SNP positions
#                 snp_positions = self.find_snps_in_consensus(ref_seq, consensus_seq)
#                 # Check if we found SNPs near expected positions (within 10bp due to deletion shifting)
#                 expected_snps_found = 0
#                 for exp_pos in expected['expected_snp_positions']:
#                     # Adjust for deletion if it's before SNP position
#                     deletion_pos = 500  # From generation script
#                     if exp_pos > deletion_pos:
#                         adjusted_pos = exp_pos - expected['expected_deletion_size']
#                     else:
#                         adjusted_pos = exp_pos
                    
#                     # Check if SNP found near expected position (within 5bp)
#                     found = any(abs(pos - adjusted_pos) <= 5 for pos in snp_positions)
#                     if found:
#                         expected_snps_found += 1
                
#                 snps_match = expected_snps_found >= 2  # At least 2 of 3 SNPs should be found
#                 validation_details.append(f"SNPs found: {len(snp_positions)} at positions {snp_positions[:5]}")
#                 validation_details.append(f"Expected SNPs found: {expected_snps_found}/3 - {'✓' if snps_match else '✗'}")
                
#                 # Test passes if length matches OR SNPs match
#                 passed = length_match or snps_match
                
#             elif test_name == "COMBINED_INDELS":
#                 # Check length difference (should be -12bp: +8 insertion - 20 deletion)
#                 length_match = abs(length_diff - expected['expected_length_diff']) <= 2  # Allow 2bp tolerance
#                 validation_details.append(f"Length difference: {length_diff} bp (expected: {expected['expected_length_diff']} bp) - {'✓' if length_match else '✗'}")
                
#                 # For INDEL-only test, only check length
#                 passed = length_match
            
#             # Write metadata
#             metadata_file = output_directory / f"{test_name}_metadata.txt"
#             with open(metadata_file, 'w') as f:
#                 f.write(f"Test: {test_name}\n")
#                 f.write(f"Description: {test_case['description']}\n")
#                 f.write(f"Expected: {test_case['expected']}\n")
#                 f.write(f"Reference length: {REF_LENGTH} bp\n")
#                 f.write(f"Consensus length: {consensus_length} bp\n")
#                 f.write(f"Length difference: {length_diff:+d} bp\n")
#                 f.write(f"\nValidation:\n")
#                 for detail in validation_details:
#                     f.write(f"  {detail}\n")
#                 f.write(f"\nResult: {'✓ PASS' if passed else '✗ FAIL'}\n")
#                 f.write(f"\nConsensus file: {output_file.name}\n")
            
#             print(f"\n{test_name}:")
#             print(f"  Consensus: {consensus_length} bp ({length_diff:+d} bp from reference)")
#             for detail in validation_details:
#                 print(f"  {detail}")
#             print(f"  Result: {'✓ PASS' if passed else '✗ FAIL'}")
            
#             # Assert test passes
#             assert passed, f"Test failed: Expected length difference {expected.get('expected_length_diff', 'N/A')} bp, got {length_diff} bp. {', '.join(validation_details)}"
            
#         finally:
#             os.chdir(original_cwd)


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])

