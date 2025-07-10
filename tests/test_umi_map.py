import unittest
from belt.scripts import umi_map
class TestExtract(unittest.TestCase):
    def test_umi_map(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250-long.bam"
        minReadLength = 2000
        motif = "NHHHHHHHHHHHHHM"
        whitelist = "../test_data/0-library_whitelist/LP250-short_correctedUMIs.csv"
        outputDir = "../test_data/5-barcode_extraction/"

        try:
            umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")

    def test_umi_map_sample2(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250v1.bam"
        minReadLength = 2000
        motif = "NHHHHHHHHHHHHHM"
        whitelist = "../test_data/0-library_whitelist/LP250-short_correctedUMIs.csv"
        outputDir = "../test_data/5-barcode_extraction/"

        try:
            umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")
