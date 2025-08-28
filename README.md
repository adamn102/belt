BELT - BalazsLab Error-correction Longread Tools

BELT is a toolkit designed for processing and analyzing long-read sequencing data. It provides functionality for extracting reads based on motifs and UMIs (Unique Molecular Identifiers), grouping reads by UMIs, and analyzing UMI distributions. BELT also allows for the integration of short-read sequencing data to error-correct UMIs from the long-read sequencing data.

## Features
- **umi-map**: Extracts UMIs from mapped reads (.bam) using a given UMI motif, corrects potential UMI errors, and returns a .csv file containing read names and their associated UMI sequences. Currently supports dually barcoded reference sets where both UMIs share the same UMI pattern (e.g. NNNNNNNNN)
- **umi** (Update to umi-group): Groups reads by their UMI sequences and a sorted UMI tag BAM files with each read tagged by a specific UMI


- **Whitelist**: Creates an error-corrected UMI database from short-read sequencing data. This database serves as a reference for correcting UMIs in long-read data.
- **Correct**: Performs error correction on long-read UMIs by comparing them against the whitelist database generated from short-read data, improving UMI accuracy.
- **Hamming**: Analyzes UMI similarity by calculating Hamming distances between UMI sequences. Includes visualization tools to display the distribution of these distances.

---

## Installation
1. Clone the repository:
```
git clone https://github.com/adamn102/belt
```

2. Install dependencies, and set up environment
```
echo "Pray to anaconda gods, May the conda be on your side, and pip away from you!"
conda env create -f [env.yaml](http://_vscodecontentref_/0)
```

If you already have an environment named belt then do this first:
If you are on the belt environment activate it first

```
conda deactivate
```

```
conda remove -n belt --all
```


### Old Notes 
TODO: sort through these later
```
python3 belt.py extract input.fastq --config config.yaml
```

Alternatively you can override anything in the config file by adding the flag as well. Eg:

```
python3 belt.py extract reads.fastq --config config.yaml --fuzzy 4
```

Note to nartsav: /Users/nomur/GitHub/test_data_belt/test_data/filtered_fastq_files/NL43-test.fastq

python3 belt.py extract /Users/nomur/GitHub/test_data_belt/test_data/raw_fastq_files/NL43-test.fastq --config config.yaml

python3 belt.py extract /Users/nartsavomur/Desktop/tmp_server_files/tests/test_data/raw_fastq_files/NL43-test.fastq --config config.yaml

ref = /Users/nomur/GitHub/test_data_belt/ERR-Seq-v2.4-NL4-3-Env-Plus-Reference.fasta

python3 belt.py extract --reference /Users/nartsavomur/Desktop/tmp_server_files/tests/ESV24-NL43-EnvPlus-Library-Quantized/ESV24-NL43-EnvPlus-Library-Reference.fasta --alignment /Users/nartsavomur/Desktop/tmp_server_files/tests/ESV24-NL43-EnvPlus-Library-Quantized/ESV24-NL43-EnvPlus-Library-sorted.bam --output /Users/nartsavomur/Desktop/tmp_server_files/tests/ESV24-NL43-EnvPlus-Library-Quantized
