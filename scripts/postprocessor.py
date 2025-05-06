import sys

import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np
from safetensors import safe_open
from scipy.stats import entropy
import torch


def get_ic(probs):
    return 2 - entropy(probs.to(torch.float32).numpy(), axis=1, base=2)


def collapse_runs(arr):
    return np.column_stack((arr[np.r_[True, np.diff(arr) > 1]], arr[np.r_[np.diff(arr) > 1, True]]))


def get_llm_exons(probs, conv_filter=[0.25, 0.5, 0.25], cutoff=0.5, min_len=10):
    llm_ic = get_ic(probs)
    ic_smooth = np.convolve(llm_ic, conv_filter, "same")
    ranges_raw = np.argwhere(ic_smooth > cutoff).flatten()
    exon_ranges = collapse_runs(ranges_raw)
    return exon_ranges[exon_ranges[:,1] - exon_ranges[:,0] > 10]


def ranges_to_coords(ranges, y):
    return [((a, y), (b, y)) for a,b in ranges]


def main() -> None:
    llm_tensors = safe_open(sys.argv[1], "pt")
    
    print("##gff-version 3")

    for i, key in enumerate(llm_tensors.keys()):
        llm_probs = llm_tensors.get_tensor(key)
        exon_ranges = get_llm_exons(llm_probs)
        for j, (start,end) in enumerate(exon_ranges):
            print(f"{key}\tllm\tCDS\t{start}\t{end}\t50\t+\t.\tID=exon_{i}_{j}")


if __name__ == "__main__":
    main()
