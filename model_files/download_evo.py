import torch
from evo import Evo

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

evo_model = Evo("evo-1.5-8k-base")
print(evo_model)
