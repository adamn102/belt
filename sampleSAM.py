import random

input_sam = "/Users/nartsavomur/Desktop/tmp_server_files/ESV24-NL43-EnvPlus-Library/ESV24-NL43-EnvPlus-Library.sam"
output_sam = "/Users/nartsavomur/Desktop/tmp_server_files/tests/ESV24-NL43-EnvPlus-Library-Quantized/ESV24-NL43-EnvPlus-Library.sam"
sample_fraction = 0.1  # 10%

with open(input_sam, 'r') as infile, open(output_sam, 'w') as outfile:
    for line in infile:
        if line.startswith('@'):
            # Copy SAM headers
            outfile.write(line)
        else:
            # Keep only ~10% of alignments
            if random.random() < sample_fraction:
                outfile.write(line)