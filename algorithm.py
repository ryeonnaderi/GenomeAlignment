from Bio import SeqIO

# Read sequences from files
lambdaRecord = SeqIO.read("Lambda.fna", "fasta")
T4phageRecord = SeqIO.read("T4phage.fna", "fasta")

# Extract sequences as strings
seq1 = str(lambdaRecord.seq)
seq2 = str(T4phageRecord.seq)

import numpy as np

def smith_waterman_gotoh(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_open_penalty=-1, gap_extension_penalty=-0.5):
    m = len(seq1)
    n = len(seq2)

    # Initialize score matrix and traceback matrix
    score_matrix = np.zeros((m+1, n+1))
    traceback_matrix = np.zeros((m+1, n+1))

    # Initialize max score and its position
    max_score = 0
    max_pos = (0, 0)

    # Fill in the score matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_open_penalty + gap_extension_penalty * traceback_matrix[i-1][j]
            insert = score_matrix[i][j-1] + gap_open_penalty + gap_extension_penalty * traceback_matrix[i][j-1]
            
            score_matrix[i][j] = max(0, match, delete, insert)
            traceback_matrix[i][j] = np.argmax([0, match, delete, insert])

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback to find the alignment
    align1 = ''
    align2 = ''
    i, j = max_pos
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == 1:  # Match or mismatch
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:  # Gap in seq1
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:  # Gap in seq2
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return max_score, align1, align2

# Example usage:

score, alignment1, alignment2 = smith_waterman_gotoh(seq1, seq2)
print("Optimal Alignment Score:", score)

