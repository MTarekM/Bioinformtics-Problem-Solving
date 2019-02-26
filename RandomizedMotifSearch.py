import random
import sys

def indexToAcid(i):
    return "ACGT"[i]

def acidToIndex(a):
    return "ACGT".index(a)


def repeatedRandomizedMotifSearch( k, t, dna ):
    bestScore = float('inf')
    bestMotifs = []
    i = 0
    while True:
        motifs = randomizedMotifSearch(k,t,dna)
        score = Score(motifs)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
            i = 0
        else:
            i += 1
        if i > 100:
            break
    return bestMotifs

def randomizedMotifSearch( k, t, dna ):
    bestMotifs = randomMotifs(dna,k)
    bestScore = Score(bestMotifs)
    while True:
        profile = profileFromMotifs(bestMotifs,1)
        motifs = motifsFromProfile(dna,profile)
        score = Score(motifs)
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
        else:
            return bestMotifs

def Score( motifs ):
    score = 0
    for count in countsFromMotifs(motifs,0):
        score += sum(count) - max(count)
    return score

def motifsFromProfile( dna, profile ):
    return [bestKmerForProfile(seq, profile) for seq in dna]

def bestKmerForProfile( seq, profile ):
    k = len(profile)
    bestProb = -1
    bestKmer = ''
    for start in range(len(seq)-k+1):
        kmer = seq[start:start+k]
        prob = probKmerInProfile(kmer, profile)
        if prob > bestProb:
            bestProb = prob
            bestKmer = kmer
    return bestKmer

def probKmerInProfile( kmer, profile ):
    prob = 1.0
    for i in range(len(kmer)):
        prob *= profile[i][acidToIndex(kmer[i])]
    return prob

def profileFromMotifs( motifs, initCount ):
    counts = countsFromMotifs(motifs,initCount)
    profile = []
    for count in counts:
        total = float(sum(count))
        probs = [c/total for c in count]
        profile.append(probs)
    return profile

def countsFromMotifs( motifs, initCount ):
    k = len(motifs[0])
    for motif in motifs:
        assert k == len(motif)
    counts = []
    for i in range(k):
        currCount = [initCount] * 4
        for motif in motifs:
            currCount[acidToIndex(motif[i])] += 1
        counts.append(currCount)
    return counts

def randomMotifs( dna, k ):
    return [randomKmer(seq,k) for seq in dna]

def randomKmer( seq, k ):
    start = random.randint(0, len(seq)-k)
    return seq[start:start+k] 





# first, import the random package

# Next, copy your RandomizedMotifSearch function (along with all required subroutines)
# from Motifs.py below this line

# Copy the ten strings occurring in the hyperlinked DosR dataset below.
dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
# set t equal to the number of strings in Dna, k equal to 15, and N equal to 100.
t=10 
k=15
N=100
# Call RandomizedMotifSearch(Dna, k, t) N times, storing the best-scoring set of motifs
# resulting from this algorithm in a variable called BestMotifs
R = 1000
BestMotifs= repeatedRandomizedMotifSearch(k,t,dna)
# Print the BestMotifs variable
print(BestMotifs)
# Print Score(BestMotifs)
print (Score(BestMotifs))


