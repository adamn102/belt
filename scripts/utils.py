from datetime import datetime
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter
import plotly.graph_objects as go

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

    def whitelist_summary_statistics(self, stats: dict):
        """
        Logs detailed UMI processing summary statistics.
        :param stats: Dictionary containing key summary statistics.
        """
        current_time_str = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
        msg_lines = [f"[BELT][{current_time_str}] Summary statistics:"]

        for key, value in stats.items():
            msg_lines.append(f"{key}: {value}")
        
        full_msg = "\n".join(msg_lines)
        print(full_msg)
        self.log_history.append(full_msg)
     
class calculator:

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
    
    def plot_hamming_distances(self, distance_dict, output_path="output_files/hamming_distances.html",sample_size=100):
        
        distances = sorted(distance_dict.keys())
        frequencies = [distance_dict[d] for d in distances]

        fig = go.Figure(
            data=[
                go.Bar(
                    x=distances,
                    y=frequencies,
                    marker_color='skyblue',
                    text=frequencies,
                    textposition='outside'
                )
            ]
        )
        fig.update_layout(
            title=f"Distribution of Hamming Distances(Sample of {sample_size})",
            xaxis_title="Hamming Distance",
            yaxis_title="Frequency",
            bargap=0.2,
            template="plotly_white"
        )
        fig.write_html(output_path)

    def plot_umi_frequencies(self, umi_list, output_path="output_files/umi_freqs.html"):
        """
        This function generates a line graph of UMI frequencies from most frequent to least frequent using Plotly.
        :param umi_list: List of UMIs (strings)
        :return: None
        """
        # Count the frequency of each UMI
        umi_counts = Counter(umi_list)
        
        # Sort UMIs by frequency in descending order
        sorted_umis = sorted(umi_counts.items(), key=lambda x: x[1], reverse=True)
        umi_labels, frequencies = zip(*sorted_umis)
        
        # Plot the line graph using Plotly
        import plotly.graph_objects as go
        fig = go.Figure(
            data=[
                go.Scatter(
                    x=list(range(len(frequencies))),
                    y=frequencies,
                    mode='lines+markers',
                    marker=dict(color='blue'),
                    line=dict(shape='linear'),
                    text=umi_labels,
                    hovertemplate='UMI: %{text}<br>Frequency: %{y}<extra></extra>'
                )
            ]
        )
        fig.update_layout(
            title='UMI Frequency Distribution',
            xaxis_title='UMIs (sorted by frequency)',
            yaxis_title='Frequency (log scale)',
            yaxis_type='log',
            template='plotly_white',
            showlegend=False
        )
        fig.write_html(output_path)

    def hamming_distance_matrix(self, umi_list, output="output_files/heatmap.html", sample_size=100):
        """
        Generates an interactive Plotly heatmap of Hamming distances between a sample of UMIs.
        The color scale is set to constant numbers (min=0, max=length of UMI).
        Hovering shows the UMIs for each cell, but axis tick labels are hidden.
        :param umi_list: List of UMIs
        :param output: Output HTML file for the heatmap
        :param sample_size: Number of UMIs to sample for the matrix
        :return: DataFrame representing the hamming distance matrix
        """
        if len(umi_list) > sample_size:
            umi_sample = np.random.choice(umi_list, size=sample_size, replace=False)
        else:
            umi_sample = umi_list

        n = len(umi_sample)
        distance_matrix = np.zeros((n, n), dtype=int)

        for i in range(n):
            for j in range(i + 1, n):
                distance = self.umi_hamming_distance(umi_sample[i], umi_sample[j])
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance

        df = pd.DataFrame(distance_matrix, index=umi_sample, columns=umi_sample)

        # Set color scale min/max based on UMI length
        umi_length = len(str(umi_sample[0])) if n > 0 else 0

        # Prepare customdata so that each cell contains [UMI_X, UMI_Y]
        customdata = np.empty((n, n, 2), dtype=object)
        for i in range(n):
            for j in range(n):
                customdata[i, j, 0] = umi_sample[j]  # x-axis UMI
                customdata[i, j, 1] = umi_sample[i]  # y-axis UMI

        fig = go.Figure(
            data=go.Heatmap(
            z=df.values,
            x=list(range(n)),  # Use integer indices for axis
            y=list(range(n)),
            colorscale='Viridis',
            colorbar=dict(title='Hamming Distance'),
            showscale=True,
            zmin=0,
            zmax=umi_length,
            customdata=customdata,
            hovertemplate="UMI X: %{customdata[0]}<br>UMI Y: %{customdata[1]}<br>Distance: %{z}<extra></extra>"
            )
        )
        fig.update_layout(
            title=f'Hamming Distance Matrix (Sample of {sample_size})',
            xaxis=dict(title='UMIs', showticklabels=False),
            yaxis=dict(title='UMIs', showticklabels=False)
        )
        fig.write_html(output)
        return df

    def chao1(self, umi_list):
        """
            Calculates the Chao1 diversity index for a list of UMIs.
            :param umi_list: List of UMIs (strings)
            :return: Chao1 diversity index (float)
            """
        counts = Counter(umi_list)
        S_obs = len(counts)
        f1 = sum(1 for c in counts.values() if c == 1)
        f2 = sum(1 for c in counts.values() if c == 2)
        if f2 == 0:
            print("Warning: No UMIs with exactly 2 counts found. Using bias-corrected formula.")
            # Avoid division by zero; use the bias-corrected formula
            return S_obs + (f1 * (f1 - 1)) / 2 if f1 > 1 else S_obs
        return S_obs + (f1 ** 2) / (2 * f2)
    
    def balazs_chao1(self, umi_list):
        """
            Calculates the updated Chao1 diversity index for a list of UMIs.
            :param umi_list: List of UMIs (strings)
            :return: Chao1 diversity index (float)
        """
        counts = Counter(umi_list)
        S_obs = len(counts)
        f1 = sum(1 for c in counts.values() if c == 1)
        f2 = sum(1 for c in counts.values() if c == 2)
        f4 = sum(1 for c in counts.values() if c == 4)
        if f4 == 0:
            print("Warning: No UMIs with exactly 4 counts found. Using bias-corrected formula.")
            # Avoid division by zero; use the bias-corrected formula
            return None 
        return S_obs + (f2 **2 ) / (2 * f4)
    
    def shannon_entropy(self, umi_list):
        """
        Calculates Shannon's entropy for a list of UMIs.
        :param umi_list: List of UMIs (strings)
        :return: Shannon's entropy (float)
        """
        counts = Counter(umi_list)
        total = sum(counts.values())
        if total == 0:
            return 0.0
        probabilities = [count / total for count in counts.values()]
        return -sum(p * np.log2(p) for p in probabilities if p > 0) 
    
    
        