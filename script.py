from random import randint

# Calculates hamming distance between two nucleotide sequences 
def hamming_dist(seqA, seqB):
    h_dist = 0 
    for i in range(0, len(seqA)):
        if seqA[i] != seqB[i]:
            h_dist += 1
    return h_dist

# Determines the consensus string from a group of motifs 
def consensus_string(motifs):
    consensus = ""
    for i in range(len(motifs[0])):
        count_A, count_C, count_G, count_T = 0, 0, 0, 0
        for motif in motifs:
            if motif[i] == "A":
                count_A += 1
            elif motif[i] == "C":
                count_C += 1
            elif motif[i] == "G":
                count_G += 1
            elif motif[i] == "T":
                count_T += 1
        if count_A >= max(count_C, count_G, count_T):
            consensus += "A"
        elif count_C >= max(count_A, count_G, count_T):
            consensus += "C"
        elif count_G >= max(count_C, count_A, count_T):
            consensus += "G"
        elif count_T >= max(count_C, count_G, count_A):
            consensus += "T"
    return consensus

# Calculates the score of a motif matrix
def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])
        score += min([hamming_dist(motif, seq*len(motif)) for seq in 'ACGT'])
    return score

# Computes the probablity of a given k-mer from a profile matrix of k-mers
def probability(text, matrix):
    prob = 1
    for i in range(0, len(text)):
        if text[i] == 'A':
            prob *= matrix['A'][i]
        elif text[i] == 'C':
            prob *= matrix['C'][i]
        elif text[i] == 'G':
            prob *= matrix['G'][i]
        else:
            prob *= matrix['T'][i]
    return prob

# Determines the profile most probable k-mer from a given DNA sequence and profile matrix 
def prof_most_prob(dna, k, prof):
    nuc_indices = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    max_prob = [-1, None]
    for i in range(len(dna)-k+1):
        current_prob = 0
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_prob *= prof[j][nuc_indices[nucleotide]]
        if current_prob > max_prob[0]:
            max_prob = [current_prob, dna[i:i+k]]
    if max_prob[1] == 0:
        return dna[0][:k]
    else:
	    return max_prob[1]

# Generates a profile matrix from a group of motifs
def generate_matrix(motifs):
    profile = {}
    A, C, G, T = [], [], [], []
    for j in range(len(motifs[0])):
        count_A, count_C, count_G, count_T = 1, 1, 1, 1
        for motif in motifs:
            if motif[j] == "A":
                count_A += 1
            elif motif[j] == "C":
                count_C += 1
            elif motif[j] == "G":
                count_G += 1
            elif motif[j] == "T":
                count_T += 1
        A.append(count_A)
        C.append(count_C)
        G.append(count_G)
        T.append(count_T)
    profile["A"] = A
    profile["C"] = C
    profile["G"] = G
    profile["T"] = T
    return profile

# Generates most probable k-mer motifs from a group of DNA sequences using a given motif profile
def generate_motifs(profile, dna, k):
    return [prof_most_prob(seq,k,profile) for seq in dna]

# Randomly generates motifs with continuously improving scores until motifs with the best score is reached
def randomized_motif_search(dna, k, t):
    random = [randint(0,len(dna[0])-k) for a in range(t)]
    motifs = [dna[i][j:j+k] for i,j in enumerate(random)]
    best_score = [score(motifs), motifs]
    while True:
        profile = generate_matrix(motifs)
        motif_group = generate_motifs(profile, dna, k)
        new_score = score(motif_group)
        if new_score < best_score[0]:
            best_score = [new_score, motif_group]
        else:
            return best_score

seqs = [
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
]

print(randomized_motif_search(seqs, 8, 5))