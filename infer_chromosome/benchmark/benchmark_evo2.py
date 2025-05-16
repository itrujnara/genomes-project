#!/usr/bin/env python3

from functools import wraps
import random
from time import perf_counter_ns

import evo2

import torch


def timeit(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        start = perf_counter_ns()
        result = f(*args, **kwargs)
        print(f"{f.__name__}: {(perf_counter_ns() - start) / 10 ** 6:.3f}ms")
        return result
    return wrapper


def make_random_seq(n, seed=42):
    nucleotides = ['A', 'C', 'G', 'T']

    random.seed(seed)

    seq = ''.join(random.choices(nucleotides, k=n))

    return seq


@timeit
def infer_no_dict(model, seq):
    with torch.inference_mode():
        return model(seq)


@timeit
def infer_with_dict(model, seq, inference_params):
    with torch.inference_mode():
        return model(seq, inference_params_dict=inference_params)


@timeit
def infer_evo_model(model, seq):
    return model(seq)


def main():
    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2("evo2_7b_base")
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)
    inference_params = model.initialize_inference_params()

    seq = make_random_seq(100)
    
    tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)
    
    print("INFERENCE MODES\n")

    logits, _ = infer_no_dict(model, tokens)

    logits, _ = infer_with_dict(model, tokens, inference_params)

    logits, _ = infer_evo_model(evo_model, tokens)

    torch.cuda.empty_cache()
    
    print()

    print("SEQUENCE LENGTHS - NO DICT\n")

    for seqlen in [100, 500, 1000, 5000, 10000]:
        print(f"Sequence length {seqlen}")
        seq = make_random_seq(seqlen)
        tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)
        logits, _ = infer_no_dict(model, tokens)
        torch.cuda.empty_cache()

    print()

    print("SEQUENCE LENGTHS - WITH DICT\n")

    for seqlen in [100, 500, 1000, 5000, 10000]:
        print(f"Sequence length {seqlen}")
        seq = make_random_seq(seqlen)
        tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)
        logits, _ = infer_with_dict(model, tokens, inference_params)
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()

