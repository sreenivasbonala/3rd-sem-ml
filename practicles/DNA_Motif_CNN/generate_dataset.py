import random
import pandas as pd

random.seed(42)

NUCLEOTIDES = "ACGT"
MOTIF = "CACGTG"
SEQ_LENGTH = 50
N_SAMPLES = 4000


def random_sequence(length):
    return "".join(random.choices(NUCLEOTIDES, k=length))


def generate_positive():
    seq = random_sequence(SEQ_LENGTH)
    position = random.randint(0, SEQ_LENGTH - len(MOTIF))
    seq = seq[:position] + MOTIF + seq[position + len(MOTIF):]
    return seq


def generate_negative():
    while True:
        seq = random_sequence(SEQ_LENGTH)
        if MOTIF not in seq:
            return seq


sequences = []
labels = []

# Positive sequences
for _ in range(N_SAMPLES // 2):
    sequences.append(generate_positive())
    labels.append(1)

# Negative sequences
for _ in range(N_SAMPLES // 2):
    sequences.append(generate_negative())
    labels.append(0)


df = pd.DataFrame({
    "sequence": sequences,
    "label": labels
})

# Shuffle dataset
df = df.sample(frac=1, random_state=42).reset_index(drop=True)

df.to_csv("dna_motif_dataset.csv", index=False)

print("====================================")
print("DNA MOTIF DATASET CREATED")
print("====================================")
print("Total sequences:", len(df))
print("Sequence length:", len(df.iloc[0]["sequence"]))
print("Motif:", MOTIF)
print("Positive samples:", (df["label"] == 1).sum())
print("Negative samples:", (df["label"] == 0).sum())

print("\nFirst 5 sequences:")
print(df.head())

print("\nDataset saved as: dna_motif_dataset.csv")
