import unittest
from belt.scripts import extract


class TestExtract(unittest.TestCase):
    def test_extract_UMIs_short(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/LP250.fa"
        alignment = "../test_data/LP250-short.bam"
        minReadLength = 250
        motif = "NHHHHHHHHHHHHHM"
        outputDir = "../test_data/test_output"

        try:
            extract.extract_UMIs(reference, alignment, minReadLength, motif, outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")

    def test_extract_UMIs_long(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/LP250.fa"
        alignment = "../test_data/LP250-long.bam"
        minReadLength = 2000
        motif = "NHHHHHHHHHHHHHM"
        outputDir = "../test_data/test_output"

        try:
            extract.extract_UMIs(reference, alignment, minReadLength, motif, outputDir)
        except Exception as e:
            self.fail(f"extract_UMIs raised an exception: {e}")
