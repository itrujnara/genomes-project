import sys

from Bio import SeqIO
from Bio.Seq import Seq

def main():
    if len(sys.argv) < 3:
        raise ValueError("Too few arguments. Usage: trim_fasta.py <in_fasta> <out_fasta>")
    
    seqs = SeqIO.parse(sys.argv[1], "fasta")
    trimmed_seqs = []

    for seq in seqs:
        seq.seq = Seq(str(seq.seq)[:8000].upper())
        trimmed_seqs.append(seq)

    SeqIO.write(trimmed_seqs, sys.argv[2], "fasta")


if __name__ == "__main__":
    main()
