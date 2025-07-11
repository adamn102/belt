import unittest
from belt.scripts import umi_map
from belt.scripts import umi_group
class TestExtract(unittest.TestCase):
    def test_umi_consensus_long(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250-long.bam"
        minReadLength = 2000
        motif = "NHHHHHHHHHHHHHM"
        grouping = "5pUMI"
        threshold = 10
        whitelist = "../test_data/0-library_whitelist/LP250-short_correctedUMIs.csv"
        outputDir = "../test_data"
        try:
            map_out_dir = outputDir + "/5-barcode_extraction"
            umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, map_out_dir)
            correctedUMIs = map_out_dir + "/LP250-long_Paired_UMIs_read_names.csv"

            group_out = outputDir + "/6-umi_group"
            umi_group.long_read_umi_grouping(reference, alignment, correctedUMIs, group_out, grouping, threshold)
        except Exception as e:
            self.fail(f"UMI consensus raised an exception: {e}")


    def test_umi_pipeline_vRNA(self):
        # Test extracting UMIs from a short read BAM file
        reference = "../test_data/0-reference_sequence/LP250.fa"
        alignment = "../test_data/4-sorted_reads/LP250v1.bam"
        minReadLength = 2000
        motif = "NHHHHHHHHHHHHHM"
        grouping = "5pUMI"
        threshold = 3
        whitelist = "../test_data/0-library_whitelist/LP250-short_correctedUMIs.csv"
        outputDir = "../test_data"
        try:
            map_out_dir = outputDir + "/5-barcode_extraction"
            umi_map.umi_extraction_map(reference, alignment, minReadLength, motif, whitelist, map_out_dir)
            correctedUMIs = map_out_dir + "/LP250v1_Paired_UMIs_read_names.csv"

            group_out = outputDir + "/6-umi_group"
            umi_group.long_read_umi_grouping(reference, alignment, correctedUMIs, group_out, grouping, threshold)
        except Exception as e:
            self.fail(f"UMI consensus raised an exception: {e}")