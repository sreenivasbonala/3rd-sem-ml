import torch
from transformers import AutoTokenizer, AutoModel

model_name = "zhihan1996/DNABERT-2-117M"

print("Loading tokenizer...")
tokenizer = AutoTokenizer.from_pretrained(
    model_name,
    trust_remote_code=True
)

print("Loading DNABERT-2 model...")
model = AutoModel.from_pretrained(
    model_name,
    trust_remote_code=True
)

print("Model loaded successfully!")

# Test DNA sequence
seq = "GTACGTACGTACGT"

# Tokenize
inputs = tokenizer(
    seq,
    return_tensors="pt"
)

print("Input shape:", inputs["input_ids"].shape)

# Generate embedding
with torch.no_grad():
    outputs = model(**inputs)

print("Output type:", type(outputs))

# DNABERT-2 returns a tuple
hidden_states = outputs[0]

print("Output shape:", hidden_states.shape)
print("Embedding dimension:", hidden_states.shape[-1])
print("Embedding generated successfully!")
