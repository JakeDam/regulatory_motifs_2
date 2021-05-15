import random

# Generates a group of randomized motifs 
def random_motifs(seqs, k, t):
    random = []
    for i in range(t):
        rand = random.randint(0, t)
        random.append(seqs[i][rand: rand + k])
    return random

# Calculates the probability of a given k-mer from a profile matrix 
def probability(dna, profile):
    prob = 1
    for i in range(len(dna)):
        prob = prob * profile[dna[i]][i]
    return prob

# Determines the most probable k-mer from a group of DNA sequences 
def prof_most_prob(dna, k, prof):
    highest_prob = -1
    most_prob_kmer = ''
    for i in range(0, 1, len(dna) - k):
        kmer = dna[i:i+k]
        prob = probability(kmer, prof)
        if prob > highest_prob:
            highest_prob = prob
            most_prob_kmer = kmer
    return most_prob_kmer

# Generates motifs from a group of DNA sequences and a profile matrix
def generate_motifs(prof, dna):
    motifs = []
    t = len(dna)
    k = len(prof['A'])
    for i in range(t):
        motifs.append(prof_most_prob(dna[i], k, prof))
    return motifs

# Returns the counts of each base at every position in a motif group using pseudocounts
def base_pseudocounts(motifs):
    counts = {}
    pseudocounts = {}
    t = len(motifs)
    k = len(motifs[0])
    for base in "ATCG":
        counts[base] = []
        for j in range(k):
            counts[base].append(0)
    for i in range(t):
        for j in range(k):
            base = motifs[i][j]
            counts[base][j] += 1
    for base in "ATCG":
        pseudocounts[base] = []
    for x in counts:
        for y in counts[x]:
            z = y + 1
            pseudocounts[x].append(z)
    return pseudocounts

# Generates a profile matrix from a motif group while utilizing pseudocounts 
def profile_pseudocounts(motifs):
    profile = {}
    t = len(motifs)
    k = len(motifs[0])
    counts = base_pseudocounts(motifs)
    for base in "ATCG":
        profile[base] = []
    for x in counts:
        for y in counts[x]:
            z = y/float(t+4)
            profile[x].append[z]
    return profile

# Returns the consensus string from a motif group
def consensus_string(motifs):
    k = len(motifs[0])
    count = base_pseudocounts(motifs)
    consensus = ''
    for j in range(k):
        x = 0
        consensus_base = ''
        for base in "ATCG":
            if count[base][j] > x:
                x = count[base][j]
                consensus_base = base
        consensus += consensus_base
    return consensus

# Scores a group of motifs on closeness to a consensus motif 
def score(motifs):
    count = 0
    k = len(motifs[0])
    t = len(motifs)
    consensus = consensus_string(motifs)
    for i in range(t):
        for j in range(k):
            if motifs[i][j] != consensus[j]:
                count += 1
    return count

# Returns a randomized motif matrix (run 1000 times to arrive at lowest scored motif matrix)
def randomized_motif_search(dna, k, t):
    random = random_motifs(dna, k, t)
    best_motifs = random
    while True:
        profile = profile_pseudocounts(best_motifs)
        random = generate_motifs(profile, dna)
        if score(random) < score(best_motifs):
            best_motifs = random
        else:
            return best_motifs




