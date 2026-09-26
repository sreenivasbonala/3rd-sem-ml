import torch
import numpy as np
import pandas as pd
from transformers import AutoTokenizer, AutoModel

# -----------------------------
# 1. Model
# -----------------------------
MODEL_NAME = "zhihan1996/DNABERT-2-117M"

print("Loading DNABERT-2 tokenizer...")
tokenizer = AutoTokenizer.from_pretrained(
    MODEL_NAME,
    trust_remote_code=True
)

print("Loading DNABERT-2 model...")
model = AutoModel.from_pretrained(
    MODEL_NAME,
    trust_remote_code=True
)

model.eval()

print("DNABERT-2 loaded successfully!")

# -----------------------------
# 2. Read FASTA
# -----------------------------
fasta_file = "ctcf_positive_clean.fasta"

sequences = []

with open(fasta_file, "r") as f:
    sequence = ""

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if sequence:
                sequences.append(sequence)
                sequence = ""
        else:
            sequence += line

    if sequence:
        sequences.append(sequence)

print("Number of sequences:", len(sequences))
print("First sequence length:", len(sequences[0]))

# -----------------------------
# 3. Generate embeddings
# -----------------------------
embeddings = []

for i, seq in enumerate(sequences):

    inputs = tokenizer(
        seq,
        return_tensors="pt"
    )

    with torch.no_grad():
        outputs = model(**inputs)

    # DNABERT-2 returns a tuple
    hidden_states = outputs[0]

    # Mean pooling across sequence tokens
    embedding = hidden_states.mean(dim=1)

    embeddings.append(
        embedding.squeeze(0).cpu().numpy()
    )

    if (i + 1) % 100 == 0:
        print(f"Processed {i + 1}/{len(sequences)} sequences")

# -----------------------------
# 4. Convert to NumPy array
# -----------------------------
embeddings = np.array(embeddings)

print("Embedding matrix shape:", embeddings.shape)

# -----------------------------
# 5. Save embeddings
# -----------------------------
np.save(
    "ctcf_dnabert2_embeddings.npy",
    embeddings
)

# Also save CSV
embedding_df = pd.DataFrame(
    embeddings,
    columns=[f"dim_{i+1}" for i in range(embeddings.shape[1])]
)

embedding_df.to_csv(
    "ctcf_dnabert2_embeddings.csv",
    index=False
)

print("Saved:")
print("  ctcf_dnabert2_embeddings.npy")
print("  ctcf_dnabert2_embeddings.csv")

print("DNABERT-2 CTCF embedding generation completed!")
