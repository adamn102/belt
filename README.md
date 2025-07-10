BELT - BalazsLab Error-correction Longread Tools

BELT is a toolkit designed for processing and analyzing long-read sequencing data. It provides functionality for extracting reads based on motifs and UMIs (Unique Molecular Identifiers), grouping reads by UMIs, and analyzing UMI distributions. BELT also allows for the integration of short-read sequencing data to error-correct UMIs from the long-read sequencing data.

## Features
- **Extract**: Extracts reads from mapped (.bam) files based on 5' and 3' motifs and UMIs and returns a .csv file containing the read name, 5' UMI sequence and 3' UMI sequence. Optionally, it will drop sequences that do not have any UMIs. As of now, this only works for dually barcoded reference sets. 
- **Group**: Groups reads by their 5' UMI sequences and outputs UMI-specific BAM files.
- **Whitelist**: Creates error corrected UMI data base from short-read seuqencing that can then be used to error correct long-read UMIs.
- **Correct**: Corrects UMIs based on short-read data, allowing for error correction of long-read UMIs.
- **Hamming**: Calculates and visualizes the distribution of Hamming distances between UMIs.

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
