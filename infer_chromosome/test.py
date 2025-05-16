#!/usr/bin/env python3

import evo2

import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

evo_model = evo2.Evo2("evo2_7b_base")

model = evo_model.model.to(device)
tokenizer = evo_model.tokenizer
inference_params = model.initialize_inference_params()

seq = "N" * 100

tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).to(device).unsqueeze(0)

logits, embeddings = evo_model(tokens, return_embeddings=True, layer_names=["norm"])

#print(logits[:,:,(65,67,71,84)])

#print(model.state_dict().keys())

print(embeddings)

print(embeddings["norm"].shape)

