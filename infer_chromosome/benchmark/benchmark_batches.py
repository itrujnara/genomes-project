#!/usr/bin/env python3

from functools import wraps
import random
from time import perf_counter_ns

from Bio import SeqIO
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
def infer_unbatched(model, seq, context_size, stride, inference_params):
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


@timeit
def infer_batched(model, seq, context_size, stride, inference_params):
    total_len = seq.shape[1]

    logits_out = torch.zeros(total_len - context_size, 512)

    batch_size = 8
    batch_subseqs = []
    batch_starts = []
    
    for start in range(0, total_len - context_size, stride):
        end = start + context_size
        subseq = seq[:, start:end]  # shape: [1, context_size]
        batch_subseqs.append(subseq)
        batch_starts.append(start)

        if len(batch_subseqs) == batch_size:
            batch = torch.cat(batch_subseqs, dim=0)  # shape: [10, context_size]
            logits, _ = infer_seq(model, batch, inference_params)  # shape: [10, context_size, 512]

            for i in range(batch_size):
                slice_start = batch_starts[i]
                slice_end = slice_start + stride
                logits_out[slice_start:slice_end, :] = logits[i, -stride:, :]

            batch_subseqs = []
            batch_starts = []

            del logits
            torch.cuda.empty_cache()

    # Handle leftovers
    if batch_subseqs:
        batch = torch.cat(batch_subseqs, dim=0)
        logits, _ = infer_seq(model, batch, inference_params)

        for i in range(len(batch_subseqs)):
            slice_start = batch_starts[i]
            slice_end = slice_start + stride
            logits_out[slice_start:slice_end, :] = logits[i, -stride:, :]

        del logits
        torch.cuda.empty_cache()

    return logits_out


def main():
    # run parameters
    context_size = 8192
    stride = 128

    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2("evo2_7b_base")
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)
    inference_params = model.initialize_inference_params()

    seqs = SeqIO.parse("/users/cn/itrujnara/genomes-project/infer_chromosome/random_chunk.fa", "fasta") # TODO make it an arg
    seq = str(list(seqs)[0].seq)

    tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)

    #logits = infer_unbatched(model, tokens, context_size, stride, inference_params)

    logits = infer_batched(model, tokens, context_size, stride, inference_params)

    tensors = {"logits": logits}

    save_file(tensors, "tensors/evo_batched.safetensors")


if __name__ == "__main__":
    main()

