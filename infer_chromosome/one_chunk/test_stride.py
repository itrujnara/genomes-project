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


def infer_seq(model, seq, inference_params):
    with torch.inference_mode():
        return model(seq, inference_params_dict=inference_params)


@timeit
def infer_sliding(model, tokenizer, seq, context_size, stride, device):
    model.eval()
    total_len = len(seq)

    logits_out = torch.zeros(total_len - context_size, 512)

    for start in range(0, total_len - context_size, stride):
        inference_params = model.initialize_inference_params()
        end = start + context_size + stride
        subseq = seq[start:end]
        tokens = torch.tensor(tokenizer.tokenize(subseq), dtype=torch.int).to(device).unsqueeze(0)

        logits, _ = infer_seq(model, tokens, inference_params)
        logits_out[start:start+stride,:] = logits[:,-stride:,:].squeeze(0)
        
        del logits, tokens

    torch.cuda.empty_cache()

    return logits_out


def main():
    torch.cuda.empty_cache()

    # run parameters
    context_size = 8192

    device = torch.device("cuda" if torch.cuda.is_available else "cpu")

    evo_model = evo2.Evo2("evo2_7b_base")
    tokenizer = evo_model.tokenizer
    model = evo_model.model.to(device)

    seqs = SeqIO.parse("/users/cn/itrujnara/genomes-project/infer_chromosome/random_chunk.fa", "fasta") # TODO make it an arg
    seq = str(list(seqs)[0].seq)

    for stride in [1,16,32,64,128,256,512,1024]:
        logits = infer_sliding(model, tokenizer, seq, context_size, stride, device)

        tensors = {"logits": logits}

        save_file(tensors, f"strides/longseq_stride{stride}.safetensors")

        del logits, tensors
        torch.cuda.empty_cache()

    print(logits.shape)
    print(logits[:10,:5])


if __name__ == "__main__":
    main()

