import unittest
from belt.scripts import consensus
class TestExtract(unittest.TestCase):
    def test_umi_consensus_long(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250-long.bam"
        correctedUMIs = "../test_data/5-barcode_extraction/LP250-long_Paired_UMIs_read_names.csv"
        outputDir = "../test_data/6-barcode_consensus/"
        grouping = "5pUMI"
        threshold = 5
        try:
            consensus.long_read_consensus(reference, alignment, correctedUMIs, outputDir, grouping,threshold)
        except Exception as e:
            self.fail(f"UMI consensus raised an exception: {e}")

    def test_umi_consensus_vRNA(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250v1.bam"
        correctedUMIs = "../test_data/5-barcode_extraction/LP250v1_Paired_UMIs_read_names.csv"
        outputDir = "../test_data/6-barcode_consensus/"
        grouping = "5pUMI"
        threshold = 5
        try:
            consensus.long_read_consensus(reference, alignment, correctedUMIs, outputDir, grouping,threshold)
        except Exception as e:
            self.fail(f"UMI consensus raised an exception: {e}")
