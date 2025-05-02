import random

class synthetic_data:
    def __init__(self, num_reads=1000, read_length=150, umi_length=12 ,
                 constant1 = "CTAGGGATAACAGGGTAATCCCAATCGTTCACCGTCGGAT",
                 constant2 = "TGTTGGCCTCGTGCATCCATGGTTGAGGTGTGAGTGA",
                 motif1 = "TCACG",
                 motif2 = "CTCGA"):
        self.num_reads = num_reads
        self.read_length = read_length
        self.umi_length = umi_length
        self.sequences = []
        self.umis = []
        self.reads = []
        self.umi_pairs = []
        self.constant1 = constant1
        self.constant2 = constant2
        self.motif1 = motif1
        self.motif2 = motif2
        
        
    def generate_umi(self):
        """
        Generate random DNA sequences of a given length.
        """
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(bases) for _ in range(self.umi_length))
    
    def generate_read(self):
        """
        Generate random DNA sequences of a given length.
        """
        umi5p = self.generate_umi()
        umi3p = self.generate_umi()
        with open("/Users/nartsavomur/Desktop/tmp_server_files/tests/nl43_env+.txt", 'r') as file:
            ref_seq = file.read()

        return self.constant1 + umi5p + self.motif1 + ref_seq + self.motif2 + self.constant2 + umi3p
    

if __name__ == "__main__":
    # Example usage
    synthetic_data_gen = synthetic_data(num_reads=1000, read_length=150, umi_length=12)
    for _ in range(10):
        print(synthetic_data_gen.generate_read())