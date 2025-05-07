BELT - BalazsLab Error-correction Longread Tools

BELT is a toolkit designed for processing and analyzing long-read sequencing data. It provides functionality for extracting reads based on motifs and UMIs (Unique Molecular Identifiers), grouping reads by UMIs, and analyzing UMI distributions.

## Features
- **Extract Reads**: Extracts reads from FASTQ files based on 5' and 3' motifs and UMIs.
- **Group Reads by UMIs**: Groups reads by their 5' UMI sequences and outputs UMI-specific BAM files.
- **Analyze UMI Hamming Distances**: Calculates and visualizes the distribution of Hamming distances between UMIs.

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

If you already have a enverioment named belt then do this first:
If you are on the belt environment deactoivate it first

```
conda deactivate
```

```
conda remove -n belt --all
```











### Old Notes 
TODO: sort through these later
'''
python3 belt.py extract input.fastq --config config.yaml
'''

Alternatively you can override anything in the config file by adding the flag as well. Eg:

'''
python3 belt.py extract reads.fastq --config config.yaml --fuzzy 4
'''

Note to nartsav: /Users/nomur/GitHub/test_data_belt/test_data/filtered_fastq_files/NL43-test.fastq

python3 belt.py extract /Users/nomur/GitHub/test_data_belt/test_data/raw_fastq_files/NL43-test.fastq --config config.yaml

python3 belt.py extract /Users/nartsavomur/Desktop/tmp_server_files/tests/test_data/raw_fastq_files/NL43-test.fastq --config config.yaml