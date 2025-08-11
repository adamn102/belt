import unittest
from belt.scripts import whitelist
from belt.scripts import whitelist_longread
class TestExtract(unittest.TestCase):
    def test_extract_UMIs(self):
        # Test extracting UMIs from a short read BAM file
        reference = "/Users/nomur/Desktop/test-data/LP250.fa"
        alignment = "/Users/nomur/Desktop/LP250.bam"
        minReadLength = 250
        #motif = "NNNNNNNNNNNN"
        motif = "NHHHHHHHHHHHHHM"
        clusterMethod = "default"
        outputDir = "/Users/nomur/Desktop/test-data/output"

        try:
            whitelist_longread.umi_whitelist(reference, alignment, minReadLength, motif, outputDir ,clusterMethod)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")
