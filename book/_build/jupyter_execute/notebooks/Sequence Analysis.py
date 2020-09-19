## Data Load / Read 

- Explore 
- Reads 
- Frequency
- Percentage 
- Use all 
- Visualization

!head ../data/den1.fasta

!head ../data/den2.fasta

!head ../data/den3.fasta

!head ../data/den4.fasta

def readFASTA(inputfile): 
    """Reads a sequence file and returns as string"""
    with open(inputfile, "r") as seqfile:
        # skip the name line 
        seq = seqfile.readline()
        seq = seqfile.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq 

# read 
seq1 = readFASTA('../data/den1.fasta')
seq2 = readFASTA('../data/den2.fasta')
seq3 = readFASTA('../data/den3.fasta')
seq4 = readFASTA('../data/den4.fasta')

## Length

print("Length of DEN1: ", len(seq1))
print("Length of DEN2: ", len(seq2))
print("Length of DEN3: ", len(seq3))
print("Length of DEN4: ", len(seq4))

## Frequency

from collections import Counter
def basecount_fast(seq): 
    """"Count the frequencies of each bases in sequence including every letter""" 
    freqs = Counter(seq)
    return freqs

print("Frequency of DEN1: ", basecount_fast(seq1))
print("Frequency of DEN2: ", basecount_fast(seq2))
print("Frequency of DEN3: ", basecount_fast(seq3))
print("Frequency of DEN4: ", basecount_fast(seq4))

def basecount(seq, useall=False, calfreqs_pc=False): 
    """Count the frequencies of each bases in sequence including every letter"""
    length = len(seq)
    
    if calfreqs_pc: 
        freq_pc = {}
    else: 
        base_counts = {} 
    
    if useall: 
        seqset = set(seq) 
    else: 
        seqset = ("A", "T", "C", "G") 
    
    for letter in seqset: 
        num = seq.count(letter)
        if calfreqs_pc:
            pc = num / length
            freq_pc[letter] = pc
        else: 
            base_counts[letter] = num 
    
    if calfreqs_pc: 
        return freq_pc
    else: 
        return base_counts

print("Frequency of DEN1: ", basecount(seq1, useall=True))
print("Frequency of DEN2: ", basecount(seq2,  useall=True))
print("Frequency of DEN3: ", basecount(seq3,  useall=True))
print("Frequency of DEN4: ", basecount(seq4,  useall=True))

print("Percentage of DEN1: ", basecount(seq1, calfreqs_pc=True))
print("Percentage of DEN2: ", basecount(seq2, calfreqs_pc=True))
print("Percentage of DEN3: ", basecount(seq3, calfreqs_pc=True))
print("Percentage of DEN4: ", basecount(seq4, calfreqs_pc=True))

f1 = basecount_fast(seq1)
f2 = basecount_fast(seq2)
f3 = basecount_fast(seq3)
f4 = basecount_fast(seq4)

import matplotlib.pyplot as plt
import seaborn as sns 
plt.rcParams['figure.figsize'] = (8,6) 
plt.rcParams['font.size'] = 14 

plt.style.use('seaborn-whitegrid')

plt.bar(f1.keys(), f1.values())
plt.title("Nucleotide Frequency Distribution of DEN1")
plt.xlabel("Bases")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show() 

plt.bar(f2.keys(), f2.values())
plt.title("Nucleotide Frequency Distribution of DEN2")
plt.xlabel("Bases")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show() 

plt.bar(f3.keys(), f3.values())
plt.title("Nucleotide Frequency Distribution of DEN3")
plt.xlabel("Bases")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show() 

plt.bar(f4.keys(), f4.values())
plt.title("Nucleotide Frequency Distribution of DEN4")
plt.xlabel("Bases")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show() 

plt.pie(f1.values(), labels=f1.keys(), autopct='%1.1f%%', shadow=True)
plt.tight_layout() 
plt.show()

plt.pie(f2.values(), labels=f2.keys(), autopct='%1.1f%%', shadow=True)
plt.tight_layout() 
plt.show()

plt.pie(f3.values(), labels=f3.keys(), autopct='%1.1f%%', shadow=True)
plt.tight_layout() 
plt.show()

plt.pie(f4.values(), labels=f4.keys(), autopct='%1.1f%%', shadow=True)
plt.tight_layout() 
plt.show()

fig, ax = plt.subplots(nrows=2, ncols=2)
print(ax) 

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True)
caption = "Figure: Nucleotide Frequency Distribution of Dengue"

# ax1 
ax1.bar(f1.keys(), f1.values())
ax1.set_title("DEN1")
ax1.set_ylabel("Frequency")

# ax2 
ax2.bar(f2.keys(), f2.values())
ax2.set_title("DEN2")
ax2.set_ylabel("Frequency")

# ax3
ax3.bar(f3.keys(), f3.values())
ax3.set_title("DEN3")
ax3.set_xlabel("Bases")
ax3.set_ylabel("Frequency")


# ax4 
ax4.bar(f4.keys(), f4.values())
ax4.set_title("DEN4")
ax4.set_xlabel("Bases")
ax4.set_ylabel("Frequency")

# Caption
plt.figtext(0.5, 0.000001, caption, wrap=True, horizontalalignment='center', fontsize=18)

# layout
plt.tight_layout()
# plt.savefig('../output_figs/den_plot.png')
plt.show() 

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

# ax1 
ax1.bar(f1.keys(), f1.values())
ax1.set_title("DEN1")
ax1.set_ylabel("Frequency")

# ax2 
ax2.pie(f1.values(), labels=f1.keys(), autopct='%1.1f%%', shadow=True)
ax2.set_title("DEN1")

# ax3
ax3.bar(f2.keys(), f2.values())
ax3.set_title("DEN2")
ax3.set_ylabel("Frequency")

# ax4 
ax4.pie(f2.values(), labels=f2.keys(), autopct='%1.1f%%', shadow=True)
ax4.set_title("DEN2")


# layout
plt.tight_layout()
# plt.savefig('den.pdf')
plt.show() 

def readFASTA(inputfile): 
    """Reads a sequence file and returns as string"""
    with open(inputfile, "r") as seqfile:
        # skip the name line 
        seq = seqfile.readline()
        seq = seqfile.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq 

# read 
seq1 = readFASTA('../data/den1.fasta')
seq2 = readFASTA('../data/den2.fasta')
seq3 = readFASTA('../data/den3.fasta')
seq4 = readFASTA('../data/den4.fasta')

len(seq1)

# den1 
round((seq1.count("G")+seq1.count("C")) / len(seq1) * 100, 2) 

len(seq2)

# den2
round((seq2.count("G")+seq2.count("C")) / len(seq2) * 100, 2) 

# den3
round((seq3.count("G")+seq3.count("C")) / len(seq3) * 100, 2) 

round((seq4.count("G")+seq4.count("C")) / len(seq4) * 100, 2) 

def calculateGC(seq): 
    """Take a sequence as input and calculate the GC %"""
    gc = round((seq.count("G")+seq.count("C")) / len(seq) * 100, 2) 
    return gc   

calculateGC(seq1)

def calculateAT(seq): 
    """Take a sequence as input and calculate the AT %"""
    at = round((seq.count("A")+seq.count("T")) / len(seq) * 100, 2) 
    return at   

calculateAT(seq1)

# Sliding Window Analysis 
len(seq1) / 10 

calculateGC(seq1[0:2001]) 

calculateGC(seq1[2001:4001]) 

len(seq1)

sub = len(seq1) - 2000 + 1 





range(len(seq1))

seq1[:10]

range(start=, stop=, step)


range(n) # n-1 
range(n+1) # n 

list(range(0, 10, 2 )) 

seq1[0]

list(range(0, 10, 2))

list(range(0, 10735+1, 2000))

list(range(0, len(seq1)+1, 2000)) 







def subSeqGC(seq, windowsize=2000): 
    """Returns sub-sequence GC Ration"""
    res = []
    for i in range(0, len(seq)-windowsize+1, windowsize): 
        subSeq = seq[i:i+windowsize] 
        gc = calculateGC(subSeq) 
        res.append(gc) 
    return res

subSeqGC(seq1)

len(subSeqGC(seq1)) 

subSeqGC(seq1, 500)

len(subSeqGC(seq1, 500)) 

import matplotlib.pyplot as plt
import seaborn as sns 
plt.rcParams['figure.figsize'] = (8,6) 
plt.rcParams['font.size'] = 14 
plt.style.use('seaborn-whitegrid')

gc = subSeqGC(seq1, 300)

range(len(gc))

gc = subSeqGC(seq1, 300)
plt.plot(range(len(gc)), gc)
plt.title("Sub-sequence GC Distribution")
plt.xlabel("Base-pair Position")
plt.ylabel("% GC")
plt.tight_layout() 
plt.savefig("../output_figs/sub_gc.pdf")
plt.show() 

gc1 = subSeqGC(seq1, 300)
gc1

gc2 = subSeqGC(seq2, 300)
gc2

gc1 = subSeqGC(seq1, 300)
plt.plot(range(len(gc1)), gc1, 'ro', label="DEN1")


gc2 = subSeqGC(seq2, 300)
plt.plot(range(len(gc2)), gc2, 'bs', label="DEN2")

gc3 = subSeqGC(seq3, 300)
plt.plot(range(len(gc3)), gc3, color='purple', label="DEN3", linestyle='--')


gc4 = subSeqGC(seq4, 300)
plt.plot(range(len(gc4)), gc4, color='cyan', label="DEN4")

plt.title("Sub-sequence GC Distribution")
plt.xlabel("Base-pair Position")
plt.ylabel("% GC")
plt.tight_layout() 
plt.legend() 
# plt.savefig("../output_figs/sub_gc.pdf")
plt.show() 

gc3 = subSeqGC(seq3, 300)


gc4 = subSeqGC(seq4, 300)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True)
# ax1 
gc1 = subSeqGC(seq1, 300)
ax1.plot(range(len(gc1)), gc1)
ax1.set_title("DEN1")
ax1.set_ylabel("% GC")

# ax2 
gc2 = subSeqGC(seq2, 300)
ax2.plot(range(len(gc2)), gc2)
ax2.set_title("DEN2")
ax2.set_ylabel("% GC")

# ax3
gc3 = subSeqGC(seq3, 300)
ax3.plot(range(len(gc3)), gc3)
ax3.set_title("DEN3")
ax3.set_xlabel("Base-pair Position")
ax3.set_ylabel("% GC")

# ax4 
gc4 = subSeqGC(seq4, 300)
ax4.plot(range(len(gc4)), gc4)
ax4.set_title("DEN4")
ax4.set_xlabel("Base-pair Position")
ax4.set_ylabel("% GC")

# layout
plt.tight_layout()
plt.savefig('../output_figs/den_plot.png')
plt.show() 

def subSeqAT(seq, windowsize=2000): 
    """Returns sub-sequence GC Ration"""
    res = []
    for i in range(0, len(seq)-windowsize+1, windowsize): 
        subSeq = seq[i:i+windowsize] 
        at = calculateAT(subSeq) 
        res.append(at)  
    return res

len(subSeqAT(seq1)) 

at = subSeqAT(seq1, 300)
plt.plot(range(len(at)), at)
plt.title("Sub-sequence AT Distribution")
plt.xlabel("Base-pair Position")
plt.ylabel("% AT")
plt.tight_layout() 
# plt.savefig("../output_figs/sub_gc.pdf")
plt.show() 

at1 = subSeqAT(seq1, 300)
plt.plot(range(len(at1)), at1, 'ro', label="DEN1")


at2 = subSeqAT(seq2, 300)
plt.plot(range(len(at2)), at2, 'bs', label="DEN2")

at3 = subSeqAT(seq3, 300)
plt.plot(range(len(at3)), at3, color='purple', label="DEN3", linestyle='--')


at4 = subSeqAT(seq4, 300)
plt.plot(range(len(at4)), at4, color='cyan', label="DEN4")

plt.title("Sub-sequence AT Distribution")
plt.xlabel("Base-pair Position")
plt.ylabel("% AT")
plt.tight_layout() 
plt.legend() 
plt.savefig("../output_figs/at_trend.pdf")
plt.show() 

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True)
# ax1 
at1 = subSeqAT(seq1, 300)
ax1.plot(range(len(at1)), at1)
ax1.set_title("DEN1")
ax1.set_ylabel("% AT")

# ax2 
at2 = subSeqAT(seq2, 300)
ax2.plot(range(len(at2)), at2)
ax2.set_title("DEN2")
ax2.set_ylabel("% AT")

# ax3
at3 = subSeqAT(seq3, 300)
ax3.plot(range(len(at3)), at3)
ax3.set_title("DEN3")
ax3.set_xlabel("Base-pair Position")
ax3.set_ylabel("% AT")

# ax4 
at4 = subSeqAT(seq4, 300)
ax4.plot(range(len(at4)), at4)
ax4.set_title("DEN4")
ax4.set_xlabel("Base-pair Position")
ax4.set_ylabel("% AT")

# layout
plt.tight_layout()
plt.savefig('../output_figs/den_plot.pdf')
plt.show() 

## K-mer Analysis 

A, T, C, G = Monomer 
AA, TT, CC, CG, CT = Dimer  
AAA, TTT, CTG = Codon=> Trimer 
K-mer 
k = 2 > AA, AT, TT
k = 3 > AAA 
K = 4 > AAAA 

range(start+1, stop, step)

list(range(10+1))

def buildKmers(seq, ksize): 
    # Store kmers 
    kmers = [] 
    # 
    n_mers = len(seq) - ksize + 1 
    for i in range(n_mers): 
        kmer = seq[i:i+ksize]
        kmers.append(kmer)
    return kmers   

seq = 'ATCGGTTAGGC'

len(seq)

len(seq) - 3 +1 

# N errors  

buildKmers('ATCGGTTAGGC',3) 

def readFASTA(inputfile): 
    """Reads a sequence file and returns as string"""
    with open(inputfile, "r") as seqfile:
        # skip the name line 
        seq = seqfile.readline()
        seq = seqfile.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("\t", "") 
    return seq 

seq1 = readFASTA('../data/den1.fasta')
seq2 = readFASTA('../data/den2.fasta')
seq3 = readFASTA('../data/den3.fasta')
seq4 = readFASTA('../data/den4.fasta') 

km1 = buildKmers(seq1, 2) 
km2 = buildKmers(seq2, 2) 
km3 = buildKmers(seq3, 2)
km4 = buildKmers(seq1, 2) 

from collections import Counter 
def kmerFrequency(kmers): 
    freq = Counter(kmers)
    return freq 

f1 = kmerFrequency(km1)

f1 

import matplotlib.pyplot as plt 

plt.bar(f1.keys(), f.values())
plt.show() 

f2 = kmerFrequency(km2)

f2 

plt.bar(f2.keys(), f2.values())
plt.show() 