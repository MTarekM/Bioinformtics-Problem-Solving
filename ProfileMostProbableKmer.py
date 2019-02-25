text = 'GCTCCGCAAAAGTAGGCGCTTGACGCTGCAGACGCCCCGATTGGAGGGGTTTATTCAGCCATTCATGAGGTAGAGTGGTATATAATCCCTCTACAATGCCGTATATGACCCGGCGAGCTCAAATTAATTCAGGGCATGTCATCGAAGCTAGCATCGGCCTCACCGCATCTCAGACGACTGATACCCGGATCCCAGACTCA'

k=8
profile = {
    'A':[0.36, 0.12, 0.24, 0.44, 0.24, 0.36, 0.44, 0.32],
'C':[0.28, 0.32, 0.32, 0.24, 0.36, 0.16, 0.2, 0.24],
'G':[0.2, 0.24, 0.16, 0.28, 0.12, 0.2, 0.16, 0.24],
'T':[0.16, 0.32, 0.28, 0.04, 0.28, 0.28, 0.2, 0.2]

}

def ProfileMostProbableKmer(text, k, profile):
    n = len(text)
    pr = {}
    most_prob_kmer = []
    for i in range(n-k+1):
        k_mer = text[i:i+k]
        probability = Pr(k_mer, profile)
        pr[k_mer] = probability
    m = max(pr.values())
    for key, value in pr.items():
        if pr[key] == m:
            most_prob_kmer.append(key)
    return most_prob_kmer[0]

def pr(text, profile):
    p = 1
    for i in range(len(text)):
        p = p * profile[text[i]][i]
    return p


print(ProfileMostProbableKmer(text, k, profile))
