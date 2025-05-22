#!/usr/bin/env python3

import argparse
from copy import deepcopy
from functools import wraps
import random
from time import perf_counter_ns

from Bio import SeqIO
import evo2
import polars as pl
from safetensors.torch import save_file

import torch

def make_parser():
    parser = argparse.ArgumentParser(prog="test_stride.py")

    parser.add_argument("--parquet", help="Path to the Parquet file with the sequence")
    parser.add_argument("--id", help="ID of the sequence to analyze")
    parser.add_argument("--context", type=int, default=8192, help="Context size")

    return parser


def timeit(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        start = perf_counter_ns()
        result = f(*args, **kwargs)
        print(f"{f.__name__}: {(perf_counter_ns() - start) / 10 ** 6:.3f}ms")
        return result
    return wrapper


def infer_seq(model, seq):
    with torch.inference_mode():
        return model(seq)


def main():
    # argument parsing
    parser = make_parser()
    args = parser.parse_args()
    # run parameters
    target_id = args.id
    context_size = args.context

    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2("evo2_7b_base")
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)

    df = pl.scan_parquet(args.parquet)
    df_filtered = df.filter(pl.col("id") == target_id)
    seq = df_filtered.collect()[0]["sequence"].item()[-40000-context_size:]
    tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)

    (logits, _), _ = evo_model(tokens)

    probs = logits[:,:,(65,67,71,84)].softmax(dim=-1).squeeze(0)

    print(probs[-1,:].squeeze(0))


if __name__ == "__main__":
    main()

