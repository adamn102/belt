def checkMotif(sequence, side):
    if side == "left":
        motif_index = sequence.find("AAGGTTAA")
        if motif_index != -1:
            remaining_seq = sequence[motif_index + 8:]
            if len(remaining_seq) >= 20:
                return remaining_seq[:20]
            else:
                return remaining_seq + 'N' * (20 - len(remaining_seq))
        else:
            return None
    elif side == "right":
        motif_index = sequence.rfind("AAGGTTAA")
        if motif_index != -1:
            remaining_seq = sequence[:motif_index]
            if len(remaining_seq) >= 20:
                return remaining_seq[-20:]
            else:
                return 'N' * (20 - len(remaining_seq)) + remaining_seq
        else:
            return None


def main():
    # Test cases
    test_sequences = [
        ("ACGTAAGGTTAAGCTAGCTAGTCACTAGGCTAG", "left"),  # Motif present, enough sequence after
        ("GCTACGTAAGGTTAAGCTTT", "left"),  # Motif present, not enough sequence after
        ("ACGTACGTACGT", "left"),  # No motif
        ("AATGGTCTCATGGTACATAAAGGTTAAGCTA", "right"),  # Motif present, enough sequence before
        ("AATTCCTTAAGGTTAAGCTAGCT", "right"),  # Motif present, not enough sequence before
        ("AAGGTTAAAAGCAGTGACAGGTGATGAA", "left"),  # No motif
        ("AAGGTTAAAAAGCACAGTCCGTAGGTTTAGTCAAATCAACCACATCCAGATAGCTTCGTAGTCCTACGCAGGTGATGATTGTGTGGCAATACAGCAGCTACCAATGCTGCTTGTGCTAGAAGCACAAGAGGAGGAAGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTAAGACCAATGCTCGAAAACTTGATGTGTGTCGCCTTCGTGCATCCATACTATTGAATACACAATGAACCGAAGGTTAAAAGGTTAAAAGCAGTGACAGGTGATGAA", "left"),  # No motif
    ]

    for sequence, side in test_sequences:
        result = checkMotif(sequence, side)
        print(f"\nTesting {side} side with sequence: {sequence}")
        print(f"Result: {result}")
        print("--------------------------------")

    return 0


if __name__ == "__main__":
    exit(main())


