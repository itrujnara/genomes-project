#!/usr/bin/env python3

import sys

import torch
from evo2 import Evo2

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_name = sys.argv[1] if len(sys.argv) > 1 else "evo2_1b_base"
    evo_model = Evo2(model_name)
    print(f"Loaded checkpoint {model_name} successfully.")


if __name__ == "__main__":
    main()
