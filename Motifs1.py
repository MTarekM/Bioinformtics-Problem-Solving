def Count(Motifs):
    count = {}
    k = len(Motifs[0])

    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count
A = 'A'
C = 'C'
T = 'T'
G = 'G'
Motifs = [[A,A,C,G,T,A], 
          [C,C,C,G,T,T], 
          [C,A,C,C,T,T], 
          [G,G,A,T,T,A],
          [T,T,C,C,G,G],
         ]
def Profile(Motifs):
    profile = {}
    t = len(Motifs)
    k = len(Motifs[0])
    CountMotifs = Count(Motifs)

    for symbol in "ACGT":
        profile[symbol] = []

    for x in CountMotifs:
        for y in CountMotifs[x]:
            z = y/float(t)
            profile[x].append(z)

    return profile
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
def Score(Motifs):
    count = 0
    k = len(Motifs[0])
    t = len(Motifs)
    ConsensusMotif = Consensus(Motifs)
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] != ConsensusMotif[j]:
                count += 1
    return count
print(Count(Motifs))
print(Profile(Motifs))
print(Consensus(Motifs))
print(Score(Motifs))
