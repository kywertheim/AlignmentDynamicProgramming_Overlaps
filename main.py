# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 17:19:31 2020

@author: Kenneth
"""

"""
This function reads a FASTA file.
"""
def readGN(fname):
    file = open(fname, 'r')
    seq = ''
    for line in file:
        if line[0] != '>':
            dummy = line.rstrip()
            seq += dummy
    return seq

"""
This function reads a FASTQ file.
"""
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

"""
This function calculates the edit distance between a pattern x and a text y.
Note that x is shorter than y.
"""
def editDistance(x, y):
    # Create distance matrix and fill it with zeros.
    # x denotes the pattern (rows).
    # y denotes the text (columns).
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first column of D
    for i in range(len(x)+1):
        D[i][0] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return min(D[-1])

"""
This function finds the length of the longest suffix of a which overlaps with a prefix of b.
"""
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

"""
This function first finds the k-mers of each read from a FASTQ file.
In doing so, it associates each k-mer with a set of reads containing the k-mer.
Then, for each read, it locates the set of reads associated with its suffix.
Finally, it applies overlap() to each read and every read in its associated set, ignoring comparison with itself.
"""
def overlap_fast(reads, min_length):
    kmer_sets = {}
    lengths = {}
    nodes_outgoing = []
    for read in reads:
        kmers = []
        for i in range(len(read) - min_length + 1):
            kmers.append(read[i:i+min_length])
        for kmer in kmers:
            if kmer not in kmer_sets:
                kmer_sets[kmer] = set()
            kmer_sets[kmer].add(read)
    for read in reads:
        read_suffix = read[-min_length:]
        matches = kmer_sets[read_suffix]
        outgoing = 0
        for match in matches:
            if match != read:
                length = overlap(read, match, min_length)
                if length != 0:
                    lengths[(read, match)] = length
                    outgoing = 1
        if outgoing == 1:
            nodes_outgoing.append(read)
    return lengths, nodes_outgoing
        
"""
Load the FASTA file, an excerpt of human chromosome 1.
"""            
seq = readGN('chr1.GRCh38.excerpt.fasta')

"""
Load the FASTQ file containing read sequences from Phi-X.
""" 
seq1, qual1 = readFastq('ERR266411_1.for_asm.fastq')

"""
Test editDistance(x, y).
"""
t = 'TATTGGCTATACGGTT'
p = 'GCGTATGC'
editDist = editDistance(p, t)
print(editDist)

"""
Question 1.
Find the edit distance between p and the genome.
"""
t = seq
p = 'GCTGATCGATCGTACG'
editDist = editDistance(p, t)
print(editDist)

"""
Question 2.
Find the edit distance between p and the genome.
"""
t = seq
p = 'GATTTACCAGATTGAG'
editDist = editDistance(p, t)
print(editDist)

"""
Question 3.
Find the number of overlaps (minimum length = 30) in the FASTQ file.
"""
lengths, nodes_outgoing = overlap_fast(seq1, 30)
N_pairs = len(lengths.keys())
print(N_pairs)

"""
Question 4.
Find the number of reads with a suffix that overlaps (minimum length = 30) with at least one other read.
"""
N_outgoing = len(nodes_outgoing)
print(N_outgoing)