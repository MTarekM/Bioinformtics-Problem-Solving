def freqw( str, int ):
    numb = {}
    for i in range(len(str)-int+1):
        if str[i:i+int] not in numb:
            numb[str[i:i+int]] = 0
        numb[str[i:i+int]] += 1
    maxim = max(numb.values())
    res = []
    for kmers in numb:
        if numb[kmers] == maxim:
            res.append(kmers)
    return res

str = 'AATGATCCGTTATGCCCGTGCCGCACAGAATGATCCGGTATTTTGAATGATCCGTTATGCCCGTGTATTTTGGTATTTTGGCCGCACAGGCCGCACAGAAATGAGTAGAATGATCCGTTATGCCCGTAATGATCCGGTATTTTGTTATGCCCGTAAATGAGTAGGCCGCACAGGCCGCACAGGTATTTTGGTATTTTGAAATGAGTAGGCCGCACAGGTATTTTGGCCGCACAGTTATGCCCGTTTATGCCCGTGTATTTTGAATGATCCGAAATGAGTAGGTATTTTGAATGATCCGTTATGCCCGTTTATGCCCGTGTATTTTGAAATGAGTAGGCCGCACAGAAATGAGTAGTTATGCCCGTTTATGCCCGTAATGATCCGAAATGAGTAGAATGATCCGAATGATCCGTTATGCCCGTGTATTTTGGCCGCACAGGCCGCACAGAATGATCCGAAATGAGTAGGCCGCACAGGTATTTTGGTATTTTGGCCGCACAGGTATTTTGTTATGCCCGTGTATTTTGAAATGAGTAGGTATTTTGAATGATCCGAAATGAGTAGGTATTTTGAAATGAGTAGAATGATCCGTTATGCCCGTTTATGCCCGTTTATGCCCGTGTATTTTGGCCGCACAGAAATGAGTAGGCCGCACAGGTATTTTGGTATTTTGAAATGAGTAGAATGATCCGGTATTTTGAATGATCCGTTATGCCCGTGCCGCACAGAATGATCCGGCCGCACAGGCCGCACAGGCCGCACAGGCCGCACAGGCCGCACAGAATGATCCGGCCGCACAGGCCGCACAGAATGATCCGAAATGAGTAGTTATGCCCGTAATGATCCGAAATGAGTAGGCCGCACAGGTATTTTGTTATGCCCGTGTATTTTGTTATGCCCGTAAATGAGTAGTTATGCCCGTTTATGCCCGTT'
int = 12
print (' '.join(freqw(str,int)))
