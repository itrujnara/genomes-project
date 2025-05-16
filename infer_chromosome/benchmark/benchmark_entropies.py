import evo2
from evo2_scoring import positional_entropies

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


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    seqs = [make_random_seq(10000, seq=42+i) for i in range(10)]

    evo2_model = evo2.Evo2("evo2_7b_base")
    model = evo2.model
    tokenizer = evo2.tokenizer

    entropies = positional_entropies

