import torch

def get_internal_acgt_probs(model, tokenizer, seq, device):
    """Get per-position nucleotide probabilities within a sequence. They go from position 1 to k+1."""
    tokens = torch.tensor(tokenizer.tokenize(seq), dtype=torch.int).unsqueeze(0).to(device)
    with torch.inference_mode():
        logits, _ = model(tokens)
    acgt_logits = logits[0,:,(65,67,71,84)]
    acgt_probs = torch.softmax(acgt_logits, dim=-1)
    return acgt_probs
