from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    seqs = SeqIO.parse("ch38_data/ch38_chr21.fna", "fasta")
    chr21 = list(seqs)[0]
    seq = chr21.seq

    size = 1048576
    stride = size - 8192
    starts = list(range(0,len(seq),stride))

    chunks = [seq[start:start+size] for start in starts]

    records = [SeqRecord(id=f"chunk{i}",description=f"chr21:{starts[i]}:{starts[i]+size}", seq=chunk) for i,chunk in enumerate(chunks)]

    SeqIO.write(records, "ch38_data/ch38_chr21_chunks_1mbp.fna", "fasta")


if __name__=="__main__":
    main()

