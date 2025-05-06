#!/usr/bin/env python3

import sys

from Bio import SeqIO
from evo import Evo
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
    model_name = "evo-1.5-8k-base"
    evo_model = Evo(model_name)
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)
    print(f"Loaded checkpoint {model_name} successfully.")

    with open(fasta_path) as f:
        seqs = SeqIO.parse(f, "fasta")
        seq = next(seqs)

    probs = get_internal_acgt_probs(model=model,
            tokenizer=tokenizer, 
            seq=str(seq.seq),
            device=device)

    tensors = {"probs": probs}

    save_file(tensors, out_path)


if __name__ == "__main__":
    main()
