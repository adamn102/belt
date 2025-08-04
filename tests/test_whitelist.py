import unittest
from belt.scripts import whitelist
class TestExtract(unittest.TestCase):
    def test_extract_UMIs(self):
        # Test extracting UMIs from a short read BAM file
        reference = "/Users/nomur/Desktop/test-data/LP250.fa"
        alignment = "/Users/nomur/Desktop/test-data/LP250.bam"
        minReadLength = 250
        #motif = "NNNNNNNNNNNN"
        motif = "NHHHHHHHHHHHHHM"
        clusterMethod = "default"
        outputDir = "/Users/nomur/Desktop/test-data/output"

        try:
            whitelist.umi_whitelist(reference, alignment, minReadLength, motif, clusterMethod,outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")
