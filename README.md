# belt
 BalazsLab Error-correction Longread Tools

##Extract Function
How to run this function

'''
python3 belt.py extract input.fastq --config config.yaml
'''

Alternatively you can override anything in the config file by adding the flag as well. Eg:

'''
python3 belt.py extract reads.fastq --config config.yaml --fuzzy 4
'''

Note to nartsav: /Users/nomur/GitHub/test_data_belt/test_data/filtered_fastq_files/NL43-test.fastq

python3 belt.py extract /Users/nomur/GitHub/test_data_belt/test_data/raw_fastq_files/NL43-test.fastq --config config.yaml