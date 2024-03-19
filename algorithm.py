from Bio import SeqIO

# Read sequences from files
lambdaRecord = SeqIO.read("Lambda.fna", "fasta")
T4phageRecord = SeqIO.read("T4phage.fna", "fasta")

# Extract sequences as strings
seq1 = str(lambdaRecord.seq)
seq2 = str(T4phageRecord.seq)

# Define match, mismatch, and gap penalties
match = 1
mismatch = -1
gap_penalty = -1

# Function to initialize the scoring matrix
def initialize_matrix(seq1, seq2):
    matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    for i in range(len(seq1) + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(len(seq2) + 1):
        matrix[0][j] = j * gap_penalty
    return matrix

# Function to perform global alignment using Needleman-Wunsch algorithm
def needleman_wunsch(seq1, seq2):
    # Initialize scoring matrix
    matrix = initialize_matrix(seq1, seq2)

    # Fill in the scoring matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # Calculate scores for three possible moves
            match_score = matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = matrix[i-1][j] + gap_penalty
            insert_score = matrix[i][j-1] + gap_penalty
            # Choose the maximum score
            matrix[i][j] = max(match_score, delete_score, insert_score)

    # Traceback to find the alignment
    align1 = ''
    align2 = ''
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and matrix[i][j] == matrix[i-1][j] + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif j > 0 and matrix[i][j] == matrix[i][j-1] + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1

    return align1, align2, matrix[len(seq1)][len(seq2)]

# Perform global alignment
alignment1, alignment2, score = needleman_wunsch(seq1, seq2)

# Print alignment and score
# print("Alignment 1:", alignment1)
# print("Alignment 2:", alignment2)
print("Alignment Score:", score)
