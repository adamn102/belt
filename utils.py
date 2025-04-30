from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter

class log:
    def __init__(self): # initializes the log, and starts the log history
        self.log_history = []
        self.log_history.append("Log initialized.")

    def log_message(self, msg):# disolay the message in the terminal and append it to the log history
        current_time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        msg = f"[BELT][{current_time_str}] {msg}"
        print(msg)
        self.log_history.append(msg)

    def extraction_log(self, input_fastq, output_fastq, written, skipped, no_umi5p, no_umi3p, no_umi) :  # gives out the information for the extraction, and appends it to the log history
        current_time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        # Format the log message
        msg = (
        f"[BELT][{current_time_str}] Extraction completed.\n"
        f"Input file: {input_fastq}\n"
        f"Output file: {output_fastq}\n"
        f"Reads written: {written}\n"
        f"Reads skipped: {skipped}\n"
        f"Reads with no 5' UMI: {no_umi5p}\n"
        f"Reads with no 3' UMI: {no_umi3p}\n"
        f"Reads with no UMIs: {no_umi}\n"
        )
        # Print the log message
        print(msg)

        # Append the log message to the log history
        self.log_history.append(msg)

class calculator():

    def __init__(self):
         self.calc_history = []

    def umi_hamming_distance(self, s1, s2):
        """
        This function calculates the hamming distance between two UMIs
        :param umi1: UMI 1
        :param umi2: UMI 2
        :return: Hamming distance
        """
        # return sum([1 for i, j in zip(umi1, umi2) if i != j])
        if len(s1) != len(s2):
            raise ValueError("Strings must be of the same length")
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

    def calculate_hamming_distances(self,umi_list, sample_size=100):
        """
        This function calculates the hamming distances between all pairs of UMIs in a list
        :param sample_size:
        :param umi_list: List of UMIs
        :return: List of hamming distances
        """

        distances = []
        umi_1 = []
        umi_2 = []
        umi_sample_list = np.random.choice(umi_list, size=sample_size, replace=False)
        for s1, s2 in combinations(umi_sample_list, 2):
            distances.append(self.umi_hamming_distance(s1, s2))
            umi_1.append(str(s1))
            umi_2.append(str(s2))
        #print(Counter(distances))
        out_df = pd.DataFrame({'UMI 1': umi_1, 'UMI 2': umi_2, 'Distance': distances})
        return Counter(distances)    
    
    def plot_hamming_distances(self,distance_dict):
        """
        This function plots the hamming distances
        :param distances: List of hamming distances
        :return: None
        """
        # Sort by distance (keys)
        distances = sorted(distance_dict.keys())
        frequencies = [distance_dict[d] for d in distances]
        print(distances)
        print(frequencies)

                # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(distances, frequencies, color='skyblue', edgecolor='black')
        ax.set_xlabel('Hamming Distance')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Hamming Distances')
        ax.set_xticks(distances)
        ax.grid(axis='y', alpha=0.7)

        # Add labels on top of bars
        ax.bar_label(bars, padding=3)

        plt.tight_layout()
        plt.show()

