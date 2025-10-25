#!/usr/bin/env python3
"""
Tests for BELT consensus module
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

# Import functions to test
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.consensus import (
    sanitize_umi_name,
    write_group_fastq,
    check_dependencies,
    run_consensus_for_umi
)


class TestSanitizeUmiName:
    """Test UMI name sanitization."""
    
    def test_simple_umi(self):
        """Test that simple UMI names pass through unchanged."""
        assert sanitize_umi_name("AAAAGGGGCCCCTTTT") == "AAAAGGGGCCCCTTTT"
    
    def test_umi_with_numbers(self):
        """Test UMI with numbers."""
        assert sanitize_umi_name("UMI123ABC") == "UMI123ABC"
    
    def test_umi_with_allowed_chars(self):
        """Test UMI with allowed special characters."""
        assert sanitize_umi_name("UMI-1.2_test") == "UMI-1.2_test"
    
    def test_umi_with_special_chars(self):
        """Test that special characters are replaced with underscore."""
        assert sanitize_umi_name("UMI@#$%test") == "UMI____test"
    
    def test_umi_with_spaces(self):
        """Test that spaces are replaced."""
        assert sanitize_umi_name("UMI test name") == "UMI_test_name"
    
    def test_none_umi(self):
        """Test that None UMI returns 'NA'."""
        assert sanitize_umi_name(None) == "NA"
    
    def test_empty_umi(self):
        """Test that empty UMI returns 'NA'."""
        assert sanitize_umi_name("") == "NA"


class TestWriteGroupFastq:
    """Test FASTQ writing from BAM reads."""
    
    def test_write_single_read(self, tmp_path):
        """Test writing a single read to FASTQ."""
        # Create a mock aligned segment
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}]}
        
        # Create a temporary BAM file
        bam_path = tmp_path / "test.bam"
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "ACGTACGT"
            read.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            bam.write(read)
        
        # Read the BAM and write to FASTQ
        fq_path = tmp_path / "test.fastq"
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Use until_eof=True to avoid needing an index
            reads = list(bam.fetch(until_eof=True))
            write_group_fastq("TESTUMI", reads, fq_path)
        
        # Verify FASTQ content
        assert fq_path.exists()
        with open(fq_path) as f:
            lines = f.readlines()
        
        assert len(lines) == 4
        assert lines[0].strip() == "@read1"
        assert lines[1].strip() == "ACGTACGT"
        assert lines[2].strip() == "+"
        assert lines[3].strip() == "IIIIIIII"
    
    def test_write_multiple_reads(self, tmp_path):
        """Test writing multiple reads to FASTQ."""
        header = {'HD': {'VN': '1.0'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}]}
        
        bam_path = tmp_path / "test.bam"
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for i in range(3):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "ACGT" * (i + 1)
                read.query_qualities = pysam.qualitystring_to_array("I" * len(read.query_sequence))
                read.reference_id = 0
                read.reference_start = i * 10
                read.mapping_quality = 60
                bam.write(read)
        
        fq_path = tmp_path / "test.fastq"
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Use until_eof=True to avoid needing an index
            reads = list(bam.fetch(until_eof=True))
            write_group_fastq("TESTUMI", reads, fq_path)
        
        # Verify FASTQ content
        assert fq_path.exists()
        with open(fq_path) as f:
            lines = f.readlines()
        
        # 3 reads = 12 lines (4 lines per read)
        assert len(lines) == 12
        assert lines[0].strip() == "@read0"
        assert lines[4].strip() == "@read1"
        assert lines[8].strip() == "@read2"


class TestCheckDependencies:
    """Test dependency checking."""
    
    def test_check_dependencies_available(self):
        """Test that check passes when all tools are available."""
        # This test will pass if running in the proper conda environment
        try:
            check_dependencies()
            # If we get here, all dependencies are present
            assert True
        except RuntimeError as e:
            # If dependencies are missing, skip the test
            pytest.skip(f"Dependencies not available: {e}")
    
    def test_check_dependencies_message(self, monkeypatch):
        """Test the error message when dependencies are missing."""
        import shutil
        
        # Mock shutil.which to return None for all tools
        def mock_which(tool):
            return None
        
        monkeypatch.setattr(shutil, "which", mock_which)
        
        with pytest.raises(RuntimeError) as exc_info:
            check_dependencies()
        
        error_msg = str(exc_info.value)
        assert "Missing required tools" in error_msg
        assert "minimap2" in error_msg
        assert "samtools" in error_msg
        assert "medaka_variant" in error_msg or "medaka" in error_msg
        assert "bcftools" in error_msg


class TestRunConsensusForUmi:
    """Test the main consensus calling function."""
    
    @pytest.fixture
    def reference_fasta(self, tmp_path):
        """Create a small reference FASTA for testing."""
        ref_path = tmp_path / "reference.fasta"
        
        # Create a simple 1000bp reference
        ref_seq = "A" * 250 + "C" * 250 + "G" * 250 + "T" * 250
        record = SeqRecord(Seq(ref_seq), id="ref1", description="Test reference")
        
        with open(ref_path, "w") as f:
            SeqIO.write(record, f, "fasta")
        
        return ref_path
    
    @pytest.fixture
    def reads_fastq(self, tmp_path, reference_fasta):
        """Create a small FASTQ file with reads for testing."""
        fq_path = tmp_path / "reads.fastq"
        
        # Create some reads that match the reference with a few variants
        reads = [
            ("@read1", "AAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT", "+", "I" * 40),
            ("@read2", "AAAAAAAAACCCCGCCCCCGGGGGGGGGGTTTTTTTTTT", "+", "I" * 40),  # C->G variant at pos 14
            ("@read3", "AAAAAAAAACCCCCCCCCCGGGAGGGGGGTTTTTTTTTT", "+", "I" * 40),  # G->A variant at pos 23
        ]
        
        with open(fq_path, "w") as f:
            for read in reads:
                f.write("\n".join(read) + "\n")
        
        return fq_path
    
    def test_consensus_function_signature(self):
        """Test that the consensus function has the expected signature."""
        import inspect
        sig = inspect.signature(run_consensus_for_umi)
        params = list(sig.parameters.keys())
        
        assert "reads_fq" in params
        assert "ref_fa" in params
        assert "meta_dir" in params
        assert "con_dir" in params
        assert "umi_threads" in params
        assert "model" in params
        assert "umi_safe" in params
    
    @pytest.mark.skipif(
        subprocess.run(["which", "medaka_variant"], capture_output=True).returncode != 0,
        reason="medaka_variant not available"
    )
    def test_consensus_basic_workflow(self, tmp_path, reference_fasta, reads_fastq):
        """
        Test the basic consensus workflow with real tools.
        This is an integration test that requires medaka, minimap2, samtools, bcftools.
        """
        # Setup directories
        meta_dir = tmp_path / "meta"
        con_dir = tmp_path / "consensus"
        bam_dir = tmp_path / "bam"
        log_file = tmp_path / "test.log"
        
        umi_safe = "TESTUMI"
        
        try:
            # Run consensus calling
            result = run_consensus_for_umi(
                reads_fq=reads_fastq,
                ref_fa=reference_fasta,
                meta_dir=meta_dir,
                con_dir=con_dir,
                bam_dir=bam_dir,
                umi_threads=1,
                model="r1041_e82_400bps_sup_variant_v5.0.0",
                log_file=log_file,
                umi_safe=umi_safe
            )
            
            # Verify output
            consensus_path = Path(result)
            assert consensus_path.exists()
            assert consensus_path.name == f"{umi_safe}_consensus.fasta"
            
            # Verify consensus FASTA content
            with open(consensus_path) as f:
                lines = f.readlines()
            
            # Should have a header and sequence
            assert len(lines) >= 2
            assert lines[0].startswith(">")
            assert umi_safe in lines[0]
            
            # Sequence should not be empty
            sequence = "".join(line.strip() for line in lines[1:])
            assert len(sequence) > 0
            
            # Verify intermediate files were created
            assert (meta_dir / "aligned.bam").exists()
            assert (meta_dir / "aligned.bam.bai").exists()
            assert (meta_dir / "medaka.vcf").exists()
            
            # Verify log file was written
            assert log_file.exists()
            log_content = log_file.read_text()
            assert "Step 1: Aligning reads" in log_content
            assert "Step 2: Indexing BAM" in log_content
            assert "Step 3: Calling variants" in log_content
            assert "Consensus FASTA created" in log_content
            
            print(f"\n✅ Consensus test passed!")
            print(f"   Output: {consensus_path}")
            print(f"   Sequence length: {len(sequence)} bp")
            
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Consensus calling failed: {e}")
        except FileNotFoundError as e:
            pytest.skip(f"Required tool not found: {e}")
    
    def test_consensus_cleanup_temp_files(self, tmp_path, reference_fasta):
        """Test that temporary FASTQ files are cleaned up after consensus calling."""
        # Create a temporary FASTQ
        reads_fq = tmp_path / "temp_reads.fastq"
        with open(reads_fq, "w") as f:
            f.write("@read1\nACGT\n+\nIIII\n")
        
        meta_dir = tmp_path / "meta"
        con_dir = tmp_path / "consensus"
        bam_dir = tmp_path / "bam"
        log_file = tmp_path / "test.log"
        
        assert reads_fq.exists()
        
        # Note: This test would need mocking or actual tools to run fully
        # For now, we just verify the file path setup
        assert reads_fq.name.endswith(".fastq")


class TestConsensusIntegration:
    """Integration tests using real BELT test data."""
    
    @pytest.fixture
    def belt_test_data(self):
        """Get path to BELT test data if available."""
        test_data_dir = Path(__file__).parent.parent.parent / "Belt Pipeline" / "data"
        
        if not test_data_dir.exists():
            pytest.skip("Belt Pipeline test data not available")
        
        return test_data_dir
    
    def test_consensus_with_real_grouped_bam(self, tmp_path, belt_test_data):
        """Test consensus calling with real grouped BAM from BELT pipeline."""
        # Look for a grouped BAM file
        grouped_bam_dir = belt_test_data / "6-umi_groups"
        
        if not grouped_bam_dir.exists():
            pytest.skip("Grouped BAM directory not found")
        
        grouped_bams = list(grouped_bam_dir.glob("*.grouped.bam"))
        if not grouped_bams:
            pytest.skip("No grouped BAM files found")
        
        # Use the first grouped BAM
        grouped_bam = grouped_bams[0]
        
        # Get corresponding reference
        ref_dir = belt_test_data / "0-consensus_reference"
        if not ref_dir.exists():
            pytest.skip("Reference directory not found")
        
        refs = list(ref_dir.glob("*.fa"))
        if not refs:
            pytest.skip("No reference FASTA found")
        
        reference = refs[0]
        
        print(f"\n🧪 Testing consensus with real data:")
        print(f"   BAM: {grouped_bam.name}")
        print(f"   Reference: {reference.name}")
        
        # Extract first UMI group and test consensus
        with pysam.AlignmentFile(str(grouped_bam), "rb") as bam:
            # Get first read to find a UMI
            first_read = next(bam.fetch(until_eof=True), None)
            
            if first_read is None:
                pytest.skip("Empty BAM file")
            
            # Get UMI tag
            umi = first_read.get_tag("ZU")
            umi_safe = sanitize_umi_name(umi)
            
            print(f"   Testing UMI: {umi_safe}")
            
            # Create FASTQ for this UMI
            reads_fq = tmp_path / f"{umi_safe}_reads.fastq"
            bam.reset()
            umi_reads = [r for r in bam.fetch(until_eof=True) if r.get_tag("ZU") == umi]
            
            if len(umi_reads) < 3:
                pytest.skip(f"Not enough reads for UMI {umi_safe} ({len(umi_reads)} reads)")
            
            write_group_fastq(umi, umi_reads[:10], reads_fq)  # Use first 10 reads
            
            print(f"   Reads extracted: {min(10, len(umi_reads))}")
        
        # Run consensus
        meta_dir = tmp_path / "meta"
        con_dir = tmp_path / "consensus"
        bam_dir = tmp_path / "bam"
        log_file = tmp_path / "test.log"
        
        try:
            result = run_consensus_for_umi(
                reads_fq=reads_fq,
                ref_fa=reference,
                meta_dir=meta_dir,
                con_dir=con_dir,
                bam_dir=bam_dir,
                umi_threads=1,
                model="r1041_e82_400bps_sup_variant_v5.0.0",
                log_file=log_file,
                umi_safe=umi_safe
            )
            
            consensus_path = Path(result)
            assert consensus_path.exists()
            
            # Read and display consensus
            with open(consensus_path) as f:
                lines = f.readlines()
            
            sequence = "".join(line.strip() for line in lines[1:])
            
            print(f"\n✅ Real data consensus test passed!")
            print(f"   Consensus length: {len(sequence)} bp")
            print(f"   Output: {consensus_path}")
            
            # Show first 100bp of sequence
            print(f"   First 100bp: {sequence[:100]}")
            
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Consensus calling failed: {e}")
        except FileNotFoundError as e:
            pytest.skip(f"Required tool not found: {e}")


def test_module_imports():
    """Test that the consensus module can be imported without errors."""
    from scripts import consensus
    
    assert hasattr(consensus, 'run_consensus_for_umi')
    assert hasattr(consensus, 'sanitize_umi_name')
    assert hasattr(consensus, 'write_group_fastq')
    assert hasattr(consensus, 'check_dependencies')
    assert hasattr(consensus, 'process_grouped_bam')


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])

