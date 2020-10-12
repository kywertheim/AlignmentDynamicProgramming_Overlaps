# AlignmentDynamicProgramming_Overlaps
Context: By using and modifying Python programs provided in the Coursera course entitled 'Algorithms for DNA Sequencing', I aligned DNA sequencing reads with an excerpt of human chromosome 1 by dynamic programming, finding the minimum edit distance in each case. I also found the overlaps between reads from Phi-X.

About:
1. editDistance(x, y) calculates the edit distance between a pattern x and a text y, where x is shorter than y. I used it to align DNA sequencing reads with an excerpt of human chromosome 1 by dynamic programming, finding the minimum edit distance in each case.
2. overlap(a, b, min_length) finds the length of the longest suffix of a which overlaps with a prefix of b.
3. overlap_fast(reads, min_length) first finds the k-mers of each read from a FASTQ file (reads). In doing so, it associates each k-mer with a set of reads containing the k-mer. Then, for each read, it locates the set of reads associated with its suffix. Finally, it applies overlap() to each read and every read in its associated set, ignoring comparison with itself.
4. Note that reverse complements are ignored in all cases.

Files:
1. main.py, chr1.GRCh38.excerpt.fasta, and ERR266411_1.for_asm.fastq must be in the same directory.
2. main.py should be implemented in Python 3.7.
3. chr1.GRCh38.excerpt.fasta is a FASTA file containing an excerpt of human chromosome 1.
4. ERR266411_1.for_asm.fastq contains sequencing reads from Phi-X.
