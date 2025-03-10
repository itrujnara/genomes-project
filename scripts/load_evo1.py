#!/usr/bin/env python3

import sys

import torch
from evo import Evo

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_name = sys.argv[1] if len(sys.argv) > 1 else "evo-1.5-8k-base"
    evo_model = Evo(model_name)
    print(f"Loaded checkpoint {model_name} successfully.")


if __name__ == "__main__":
    main()
