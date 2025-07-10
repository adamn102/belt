import unittest
from belt.scripts import whitelist
class TestExtract(unittest.TestCase):
    def test_extract_UMIs_short(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250-short.bam"
        minReadLength = 250
        motif = "NHHHHHHHHHHHHHM"
        clusterMethod = "default"
        outputDir = "../test_data/test_output/whitelist"

        try:
            whitelist.umi_whitelist(reference, alignment, minReadLength, motif, clusterMethod,outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")
