from datetime import datetime

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


