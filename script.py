from random import randint

# Calculates hamming distance between two nucleotide sequences 
def hamming_dist(seqA, seqB):
    h_dist = 0 
    for i in range(0, len(seqA)):
        if seqA[i] != seqB[i]:
            h_dist += 1
    return h_dist

# Calculates the score of a motif matrix
def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])
        score += min([hamming_dist(motif, seq*len(motif)) for seq in 'ACGT'])
    return score


# Determines the profile most probable k-mer from a given DNA sequence and profile matrix 
def prof_most_prob(dna, k, prof):
    bases = {nucleotide:index for index,nucleotide in enumerate('ACGT')}
    max_prob = [-1, None]
    for i in range(len(dna)-k+1):
        current_prob = 1
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_prob *= prof[j][bases[nucleotide]]
        if current_prob > max_prob[0]:
            max_prob = [current_prob, dna[i:i+k]]
    return max_prob[1]

# Returns a motif matrix profile using pseudocounts
def prof_with_pseudocounts(motifs):
	prof = []
	for i in range(len(motifs[0])):
		col = ''.join([motifs[j][i] for j in range(len(motifs))])
		prof.append([float(col.count(base)+1)/float(len(col)+4) for base in 'ACGT'])
	return prof

# Generates most probable k-mer motifs from a group of DNA sequences using a given motif profile
def generate_motifs(profile, dna, k):
    return [prof_most_prob(seq,k,profile) for seq in dna]

# Randomly generates motifs with continuously improving scores until motifs with the best score is reached
def randomized_motif_search(dna, k, t):
    random = [randint(0,len(dna[0])-k) for a in range(t)]
    motifs = [dna[i][j:j+k] for i,j in enumerate(random)]
    best_score = [score(motifs), motifs]
    while True:
        profile = prof_with_pseudocounts(motifs)
        motif_group = generate_motifs(profile, dna, k)
        new_score = score(motif_group)
        if new_score < best_score[0]:
            best_score = [new_score, motif_group]
        else:
            return best_score



with open('dataset_161_5.txt') as f:
    content = f.readlines()
seqs = [x.strip() for x in content] 


i = 0

last_motifs = (randomized_motif_search(seqs, 15, 20))

while (i < 1000):
    best_motifs = (randomized_motif_search(seqs, 15, 20))
    if score(best_motifs[1]) < score(last_motifs[1]):
        last_motifs = best_motifs[:]
    i += 1

for motif in last_motifs[1]:
    print(motif)
    print('\n')
