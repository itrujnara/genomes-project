#!/usr/bin/env python3

import sys

from Bio import SeqIO
from evo2 import Evo2
from safetensors.torch import save_file
import torch

from utils import get_internal_acgt_probs 

def main():
    # parse arguments
    if len(sys.argv) < 3:
        raise ValueError("Too few arguments. Usage: evo_inference.py <fasta> <output_path>")
    fasta_path = sys.argv[1]
    out_path = sys.argv[2]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_name = "evo2_7b"
    evo_model = Evo2(model_name)
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)
    print(f"Loaded checkpoint {model_name} successfully.")

    with open(fasta_path) as f:
        seqs = SeqIO.parse(f, "fasta")

        probs = {}

        for seq in seqs:
            probs[seq.id] = get_internal_acgt_probs(model=model,
                tokenizer=tokenizer, 
                seq=str(seq.seq),
                device=device)

    save_file(probs, out_path)


if __name__ == "__main__":
    main()
