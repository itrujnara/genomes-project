#!/usr/bin/env python3

import argparse
from copy import deepcopy
from functools import wraps
import random
import sys
from time import perf_counter_ns

from Bio import SeqIO
import evo2
import polars as pl
from safetensors.torch import save_file

import torch

def make_parser():
    parser = argparse.ArgumentParser(prog="infer_chunk.py")

    parser.add_argument("--parquet", help="Path to the Parquet file with the sequence")
    parser.add_argument("--id", help="ID of the sequence to analyze")
    parser.add_argument("--model", default="evo2_7b_base", help="Evo 2 checkpoint to use")
    parser.add_argument("--context", type=int, default=8192, help="Context size")
    parser.add_argument("--stride", type=int, default=128, help="Sliding window stride")
    parser.add_argument("--prefix", default="evo", help="Prefix for file disambiguation")

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


@timeit
def infer_sliding(model, tokenizer, seq, context_size, stride, device):
    total_len = len(seq)

    logits_out = torch.zeros(total_len - context_size, 512)
    embeddings_out = torch.zeros(total_len - context_size, 4096)

    for start in range(0, total_len - context_size, stride):
        print(f"{start}/{total_len-context_size}")
        if (start // stride) % 10 == 0:
            sys.stdout.flush()

        end = start + context_size + stride
        subseq = seq[start:end]
        tokens = torch.tensor(tokenizer.tokenize(subseq), dtype=torch.int).to(device).unsqueeze(0)
        
        (logits, _), embeddings = model(tokens, return_embeddings=True, layer_names=["norm"])
        logits_out[start:start+stride,:] = logits[:,-stride:,:].squeeze(0)
        embeddings_out[start:start+stride,:] = embeddings["norm"][:,-stride:,:].squeeze(0)
        
        del logits, tokens, embeddings

    torch.cuda.empty_cache()

    return logits_out, embeddings_out


def main():
    # argument parsing
    parser = make_parser()
    args = parser.parse_args()
    # run parameters
    target_id = args.id
    context_size = args.context
    stride = args.stride
    prefix = args.prefix

    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2(args.model)
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)

    df = pl.scan_parquet(args.parquet)
    df_filtered = df.filter(pl.col("id") == target_id)
    seq = df_filtered.collect()[0]["sequence"].item()

    logits, embeddings = infer_sliding(
                evo_model, 
                tokenizer, 
                seq, 
                context_size, 
                stride,
                device
            )

    probs = logits[:,(65,67,71,84)].softmax(dim=-1)

    tensors_probs = {"probs": probs}
    tensors_emb = {"embeddings": embeddings}

    save_file(tensors_probs, f"tensors/{prefix}_probs_context{context_size}_stride{stride}.safetensors")

    save_file(tensors_emb, f"tensors/{prefix}_embeddings_context{context_size}_stride{stride}.safetensors")

    print(logits.shape)


if __name__ == "__main__":
    main()

