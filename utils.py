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

    def write_log(self, output_file):
        """
        This function writes the log history to a file
        :param output_file: Output file name
        :return: None
        """
        with open(output_file, 'w') as f:
            for line in self.log_history:
                f.write(line + '\n')
        print(f"Log written to {output_file}")
        self.log_history.append(f"Log written to {output_file}")
        # Add a timestamp to the log
        current_time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        self.log_history.append(f"[BELT][{current_time_str}] Log written to {output_file}")

     
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
    
    def plot_hamming_distances(self,distance_dict,output_path="output_files/hamming_distances.png"):
        """
        This function plots the hamming distances
        :param distances: List of hamming distances
        :return: None
        """
        # Sort by distance (keys)
        distances = sorted(distance_dict.keys())
        frequencies = [distance_dict[d] for d in distances]

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
        plt.savefig(output_path)

    def plot_umi_frequencies(self, umi_list, output_path="output_files/umi_freqs.png"):
        """
        This function generates a line graph of UMI frequencies from most frequent to least frequent.
        :param umi_list: List of UMIs (strings)
        :return: None
        """
        # Count the frequency of each UMI
        umi_counts = Counter(umi_list)
        
        # Sort UMIs by frequency in descending order
        sorted_umis = sorted(umi_counts.items(), key=lambda x: x[1], reverse=True)
        umi_labels, frequencies = zip(*sorted_umis)
        
        # Plot the line graph
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(len(frequencies)), frequencies, marker='o', linestyle='-', color='blue')
        ax.set_xlabel('UMIs (sorted by frequency)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title('UMI Frequency Distribution', fontsize=14)
        ax.set_yscale('log')
        ax.grid(axis='y', alpha=0.7)
        
        # Optional: Limit the x-axis to show only the top N UMIs
        ax.set_xlim(0, len(frequencies))
        
        plt.tight_layout()
        plt.savefig(output_path)

    def plot_unique_3p_per_5p(self, input_fastq, output_path="output_files/concordance.png"):
        """
        This function generates a bar graph showing the number of unique 3' UMIs for each 5' UMI
        by extracting UMIs directly from the read names in a FASTQ file.
        
        :param input_fastq: Path to the input FASTQ file.
        :param output_path: Path to save the output plot.
        :return: None
        """
        from Bio import SeqIO
        from collections import defaultdict
        input_fastq = input_fastq

        # Dictionary to store 5' UMIs as keys and associated 3' UMIs as values
        umi_pair_dict = defaultdict(list)

        # Parse the FASTQ file and extract UMIs from read names
        for record in SeqIO.parse(input_fastq, "fastq"):
            if "UMI5-" in record.id and "-UMI3-" in record.id:
                umi_5p = record.id.split("UMI5-")[1].split("-UMI3")[0]
                umi_3p = record.id.split("-UMI3-")[1]
                umi_pair_dict[umi_5p].append(umi_3p)

        # Calculate the percentage of unique 3' UMIs for each 5' UMI
        unique_3p_percentage = {
            umi_5p: (len(set(umi_3p_list)) / len(umi_3p_list)) * 100
            for umi_5p, umi_3p_list in umi_pair_dict.items()
        }

        # Sort the 5' UMIs by the percentage of unique 3' UMIs in descending order
        sorted_umis = sorted(unique_3p_percentage.items(), key=lambda x: x[1], reverse=True)
        umi_5p_labels, unique_percentages = zip(*sorted_umis) if sorted_umis else ([], [])

        # Plot the bar graph
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(umi_5p_labels, unique_percentages, color='skyblue', edgecolor='black')
        ax.set_xlabel("5' UMIs", fontsize=12)
        ax.set_ylabel("Percentage of Unique 3' UMIs", fontsize=12)
        ax.set_title("Percentage of Unique 3' UMIs per 5' UMI", fontsize=14)
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # Hide x-axis labels for clarity
        ax.grid(axis='y', alpha=0.7)

        plt.tight_layout()
        plt.savefig(output_path)

    def plot_most_frequent_3p_percentage_per_5p(self, input_fastq, output_path="output_files/concordance.png"):
        """
        This function generates a bar graph showing the percentage of the most frequent 3' UMI
        for each 5' UMI by extracting UMIs directly from the read names in a FASTQ file.
        
        :param input_fastq: Path to the input FASTQ file.
        :param output_path: Path to save the output plot.
        :return: None
        """
        from Bio import SeqIO
        from collections import defaultdict, Counter

        # Dictionary to store 5' UMIs as keys and associated 3' UMIs as values
        umi_pair_dict = defaultdict(list)

        # Parse the FASTQ file and extract UMIs from read names
        for record in SeqIO.parse(input_fastq, "fastq"):
            if "UMI5-" in record.id and "-UMI3-" in record.id:
                umi_5p = record.id.split("UMI5-")[1].split("-UMI3")[0]
                umi_3p = record.id.split("-UMI3-")[1]
                umi_pair_dict[umi_5p].append(umi_3p)

        # Calculate the percentage of the most frequent 3' UMI for each 5' UMI
        most_frequent_3p_percentage = {}
        for umi_5p, umi_3p_list in umi_pair_dict.items():
            umi_3p_counts = Counter(umi_3p_list)
            most_frequent_count = max(umi_3p_counts.values())
            most_frequent_3p_percentage[umi_5p] = (most_frequent_count / len(umi_3p_list)) * 100

        # Sort the 5' UMIs by the percentage of the most frequent 3' UMI in descending order
        sorted_umis = sorted(most_frequent_3p_percentage.items(), key=lambda x: x[1], reverse=True)
        umi_5p_labels, most_frequent_percentages = zip(*sorted_umis) if sorted_umis else ([], [])

        # Plot the bar graph
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(umi_5p_labels, most_frequent_percentages, color='skyblue', edgecolor='black')
        ax.set_xlabel("5' UMIs", fontsize=12)
        ax.set_ylabel("Percentage of Most Frequent 3' UMI", fontsize=12)
        ax.set_title("Percentage of Most Frequent 3' UMI per 5' UMI", fontsize=14)
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # Hide x-axis labels for clarity
        ax.grid(axis='y', alpha=0.7)

        plt.tight_layout()
        plt.savefig(output_path)

    def plot_most_frequent_3p_percentage_and_count_per_5p(self, input_fastq, output_path="output_files/concordance.png"):
        """
        This function generates a combined stacked bar and line graph:
        - Stacked bar graph: Percentages of the most frequent, second most frequent, and third most frequent 3' UMIs for each 5' UMI.
        - Line graph: Total number of 3' UMIs for each 5' UMI.
        Additionally, it returns a dictionary mapping each 5' UMI to its most frequent 3' UMI.
        
        :param input_fastq: Path to the input FASTQ file.
        :param output_path: Path to save the output plot.
        :return: Dictionary mapping 5' UMIs to their most frequent 3' UMIs.
        """
        from Bio import SeqIO
        from collections import defaultdict, Counter

        # Dictionary to store 5' UMIs as keys and associated 3' UMIs as values
        umi_pair_dict = defaultdict(list)

        # Parse the FASTQ file and extract UMIs from read names
        for record in SeqIO.parse(input_fastq, "fastq"):
            if "UMI5-" in record.id and "-UMI3-" in record.id:
                umi_5p = record.id.split("UMI5-")[1].split("-UMI3")[0]
                umi_3p = record.id.split("-UMI3-")[1]
                umi_pair_dict[umi_5p].append(umi_3p)

        # Calculate the percentages of the most frequent, second most frequent, and third most frequent 3' UMIs
        most_frequent_3p_percentages = {}
        second_most_frequent_3p_percentages = {}
        third_most_frequent_3p_percentages = {}
        total_3p_counts = {}
        most_frequent_3p_dict = {}

        for umi_5p, umi_3p_list in umi_pair_dict.items():
            umi_3p_counts = Counter(umi_3p_list)
            sorted_counts = sorted(umi_3p_counts.items(), key=lambda x: x[1], reverse=True)

            # Get the most frequent, second most frequent, and third most frequent counts
            most_frequent_count = sorted_counts[0][1] if len(sorted_counts) > 0 else 0
            second_most_frequent_count = sorted_counts[1][1] if len(sorted_counts) > 1 else 0
            third_most_frequent_count = sorted_counts[2][1] if len(sorted_counts) > 2 else 0

            # Calculate percentages
            total_count = len(umi_3p_list)
            most_frequent_3p_percentages[umi_5p] = (most_frequent_count / total_count) * 100 if total_count > 0 else 0
            second_most_frequent_3p_percentages[umi_5p] = (second_most_frequent_count / total_count) * 100 if total_count > 0 else 0
            third_most_frequent_3p_percentages[umi_5p] = (third_most_frequent_count / total_count) * 100 if total_count > 0 else 0
            total_3p_counts[umi_5p] = total_count

            # Store the most frequent 3' UMI
            most_frequent_3p_dict[umi_5p] = sorted_counts[0][0] if len(sorted_counts) > 0 else None

        # Sort the 5' UMIs by the percentage of the most frequent 3' UMI in descending order
        sorted_umis = sorted(most_frequent_3p_percentages.items(), key=lambda x: x[1], reverse=True)
        umi_5p_labels, most_frequent_percentages = zip(*sorted_umis) if sorted_umis else ([], [])
        second_most_frequent_percentages = [second_most_frequent_3p_percentages[umi_5p] for umi_5p in umi_5p_labels]
        third_most_frequent_percentages = [third_most_frequent_3p_percentages[umi_5p] for umi_5p in umi_5p_labels]
        total_counts_sorted = [total_3p_counts[umi_5p] for umi_5p in umi_5p_labels]

        # Plot the combined graph
        fig, ax1 = plt.subplots(figsize=(12, 6))

        # Stacked bar graph for the percentages of the most frequent, second most frequent, and third most frequent 3' UMIs
        x_indices = range(len(umi_5p_labels))
        ax1.bar(x_indices, most_frequent_percentages, color='skyblue', edgecolor='black', label="Most Frequent 3' UMI (%)")
        ax1.bar(x_indices, second_most_frequent_percentages, bottom=most_frequent_percentages, color='orange', edgecolor='black', label="Second Most Frequent 3' UMI (%)")
        ax1.bar(x_indices, third_most_frequent_percentages, bottom=[i + j for i, j in zip(most_frequent_percentages, second_most_frequent_percentages)], color='green', edgecolor='black', label="Third Most Frequent 3' UMI (%)")

        ax1.set_xlabel("5' UMIs", fontsize=12)
        ax1.set_ylabel("Percentage of 3' UMIs", fontsize=12, color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')
        ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # Hide x-axis labels for clarity
        ax1.grid(axis='y', alpha=0.7)

        # Line graph for the total number of 3' UMIs
        ax2 = ax1.twinx()
        ax2.plot(x_indices, total_counts_sorted, color='red', marker='o', linestyle='-', label="Total 3' UMI Count")
        ax2.set_ylabel("Total Number of 3' UMIs", fontsize=12, color='red')
        ax2.set_yscale('log')  # Set the right y-axis to log scale
        ax2.tick_params(axis='y', labelcolor='red')

        # Add legends
        ax1.legend(loc='upper left')
        ax2.legend(loc='upper right')

        # Add title
        plt.title("Most Frequent 3' UMI Percentages and Total 3' UMI Count per 5' UMI", fontsize=14)

        plt.tight_layout()
        plt.savefig(output_path)

        # Return the dictionary of most frequent 3' UMIs
        return most_frequent_3p_dict

    def export_5p_3p_UMIs(self, input_fastq, output_path="output_files/5p_3p_UMIs.csv"):
        """
        This function extracts 5' and 3' UMIs from the read names in a FASTQ file and exports them to a CSV file.
        
        :param input_fastq: Path to the input FASTQ file.
        :param output_path: Path to save the output CSV file.
        :return: None
        """
        from Bio import SeqIO
        import pandas as pd
        from collections import defaultdict

        # Dictionary to store 5' UMIs as keys and associated 3' UMIs as values
        umi_pair_dict = defaultdict(list)

        # Parse the FASTQ file and extract UMIs from read names
        for record in SeqIO.parse(input_fastq, "fastq"):
            if "UMI5-" in record.id and "-UMI3-" in record.id:
                umi_5p = record.id.split("UMI5-")[1].split("-UMI3")[0]
                umi_3p = record.id.split("-UMI3-")[1]
                umi_pair_dict[umi_5p].append(umi_3p)

        # Prepare data for the DataFrame
        data = [{'5\' UMI': umi_5p, '3\' UMIs': ', '.join(umi_3p_list)} for umi_5p, umi_3p_list in umi_pair_dict.items()]

        # Create a DataFrame
        umi_df = pd.DataFrame(data)

        # Save the DataFrame to a CSV file
        umi_df.to_csv(output_path, index=False)     

