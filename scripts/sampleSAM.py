import random
import pysam

input_sam = "/Users/nomur/Desktop/LP250.bam"
output_sam = "/Users/nomur/Desktop/test-data/LP250_sampled.bam"
sample_fraction = 0.2  # 10%
with pysam.AlignmentFile(input_sam, "rb") as infile, \
    pysam.AlignmentFile(output_sam, "wb", template=infile) as outfile:
    for read in infile.fetch(until_eof=True):
       if random.random() < sample_fraction:
          outfile.write(read)

# Index the sampled BAM file to create a .bai file
pysam.index(output_sam)
