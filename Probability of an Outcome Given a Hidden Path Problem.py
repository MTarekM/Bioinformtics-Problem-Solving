#Given: A string x, followed by the alphabet Σ from which x was constructed, followed by a hidden path π, followed by the states States and emission matrix Emission of an HMM (Σ, States, Transition, Emission).

#Return: The conditional probability Pr(x|π) that string x will be emitted by the HMM given the hidden path π.

prob={}
with open('data.txt.txt') as f:
    list = f.readlines()
    xy = list[0].strip().upper()
    ab = list[4].strip()
    __,AX,AY,AZ = list[9].split()
    __,BX,BY,BZ = list[10].split()
    prob['AX']=float(AX)
    prob['AY']=float(AY)
    prob['AZ']=float(AZ)
    prob['BX']=float(BX)
    prob['BY']=float(BY)
    prob['BZ']=float(BZ)
total=1
for i in range(len(xy)):
    total*=prob[ab[i]+xy[i]]
print (total)
