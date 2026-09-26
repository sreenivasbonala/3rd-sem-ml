import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel

MODEL_NAME = "zhihan1996/DNABERT-2-117M"
INPUT_FASTA = "ctcf_negative_clean.fasta"

OUTPUT_NPY = "ctcf_negative_dnabert2_embeddings.npy"
OUTPUT_CSV = "ctcf_negative_dnabert2_embeddings.csv"

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

# Read negative sequences
records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

print("Number of negative sequences:", len(records))
print("First sequence length:", len(records[0].seq))

embeddings = []

for i, record in enumerate(records):

    sequence = str(record.seq)

    inputs = tokenizer(
        sequence,
        return_tensors="pt"
    )

    with torch.no_grad():
        outputs = model(**inputs)

    # DNABERT-2 returns a tuple
    hidden_states = outputs[0]

    # Mean pooling across sequence tokens
    embedding = hidden_states.mean(dim=1).squeeze(0).numpy()

    embeddings.append(embedding)

    if (i + 1) % 100 == 0:
        print(f"Processed {i + 1}/{len(records)} sequences")

# Convert to NumPy array
embeddings = np.array(embeddings)

print("Embedding matrix shape:", embeddings.shape)

# Save NumPy
np.save(OUTPUT_NPY, embeddings)

# Save CSV
df = pd.DataFrame(embeddings)
df.to_csv(OUTPUT_CSV, index=False)

print("Saved:")
print(" ", OUTPUT_NPY)
print(" ", OUTPUT_CSV)

print("DNABERT-2 negative embedding generation completed!")
