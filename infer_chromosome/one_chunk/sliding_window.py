#!/usr/bin/env python3

from functools import wraps
import random
from time import perf_counter_ns

import evo2
from safetensors.torch import save_file

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


def infer_seq(model, seq, inference_params):
    with torch.inference_mode():
        return model(seq, inference_params_dict=inference_params)


@timeit
def infer_sliding(model, seq, context_size, stride, inference_params):
    total_len = seq.shape[1]

    logits_out = torch.zeros(total_len - context_size, 512)

    for start in range(0, total_len - context_size - stride, stride):
        end = start + context_size + stride
        subseq = seq[:,start:end]
        
        logits, _ = infer_seq(model, subseq, inference_params)
        logits_out[start:start+stride,:] = logits[:,-stride:,:].squeeze(0)
        
        del logits

    torch.cuda.empty_cache()

    return logits_out


def main():
    # run parameters
    total_len = 100_000
    context_size = 10_000
    stride = 100

    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2("evo2_7b_base")
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)
    inference_params = model.initialize_inference_params()

    seq = make_random_seq(total_len)

    tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)

    logits = infer_sliding(model, tokens, context_size, stride, inference_params)

    tensors = {"logits": logits}

    save_file(tensors, "sliding_window_longseq.safetensors")

    print(logits.shape)
    print(logits[:10,:5])


if __name__ == "__main__":
    main()

