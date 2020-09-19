12/33

seq1 = 'ATGGACCAGATATAGGGAGAGCCAGGTAGGACA'
seq2 = 'ATNGACCAGATATTGNGAGAGCCGGGTAGGACN'

# length 
len(seq1)

len(seq2)

print("No. of A: ", seq1.count("A"))
print("No. of T: ", seq1.count("T"))
print("No. of G: ", seq1.count("G"))
print("No. of C: ", seq1.count("C"))

print("RF. of A: ", seq1.count("A")/len(seq1))
print("RF. of T: ", seq1.count("T")/len(seq1))
print("RF. of G: ", seq1.count("G")/len(seq1))
print("RF. of C: ", seq1.count("C")/len(seq1))

print("RF. of A: ", seq2.count("A")/len(seq2))
print("RF. of T: ", seq2.count("T")/len(seq2))
print("RF. of G: ", seq2.count("G")/len(seq2))
print("RF. of C: ", seq2.count("C")/len(seq2))
print("RF. of N: ", seq2.count("N")/len(seq2))

def basecount(seq): 
    """Count the frequencies of each bases in sequence including every letter""" 
    # {A: 3 30, C: 50, G: 4} 
    freq = {}
    for base in seq: 
        if base in freq: 
            freq[base] += 1 
        else: 
            freq[base] = 1 
    return freq     

seq_freq = basecount(seq1)

%%timeit
basecount(seq1)

seq_freq

y = seq_freq.values() 

x = seq_freq.keys() 

import matplotlib.pyplot as plt 
plt.figure(figsize=(10,6))
plt.title("Nucleotide Frequencies of DNA Sequence", fontsize=14)
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.tight_layout()
plt.bar(x, y)
plt.savefig("frequency.pdf")
plt.show() 

import matplotlib.pyplot as plt 
plt.figure(figsize=(10,6))
plt.title("Nucleotide Frequencies of DNA Sequence", fontsize=14)
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.tight_layout()
plt.bar(seq_freq.keys(), seq_freq.values())
# plt.savefig("frequency.pdf")
plt.show() 

from collections import Counter 
def counts_fast(seq): 
    """Count the frequencies of each bases in sequence including every letter""" 
    freq = Counter(seq) 
    return freq

%%timeit 
counts_fast(seq1)

seq2 = 'ATNGACCAGATATTGNGAGAGCCGGGTAGGACN'

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

basecount(seq1)

seq_info = basecount(seq1, calfreqs_pc=True)
seq_info

import matplotlib.pyplot as plt 
plt.figure(figsize=(10,6))
plt.title("Nucleotide Frequencies of DNA Sequence", fontsize=14)
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Relative Frequency", fontsize=14)
plt.tight_layout()
plt.bar(seq_info.keys(), seq_info.values())
# plt.savefig("frequency.pdf")
plt.show() 

seq2

seq2 = basecount(seq2, useall=True)

seq2

import matplotlib.pyplot as plt 
plt.figure(figsize=(10,6))
plt.title("Nucleotide Frequencies of DNA Sequence", fontsize=14)
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Relative Frequency", fontsize=14)
plt.tight_layout()
plt.bar(seq2.keys(), seq2.values())
# plt.savefig("frequency.pdf")
plt.show() 

seq1 = basecount(seq1)
seq2 = basecount(seq2)

import matplotlib.pyplot as plt 
plt.figure(figsize=(10,6))


plt.subplot(211) # row = 2 , column = 1,  figure = 1 
plt.bar(seq1.keys(), seq1.values())
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Frequency", fontsize=14)


plt.subplot(212) # row = 2,  column = 1, figure = 2 
plt.bar(seq2.keys(), seq2.values())
plt.xlabel("Bases", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.tight_layout()
plt.savefig("multi_plot.pdf")
plt.show() 