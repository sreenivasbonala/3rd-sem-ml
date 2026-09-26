import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

from Bio import SeqIO


# ============================================================
# REAL CTCF CNN - MOTIF INTERPRETABILITY
# Gradient x Input Saliency
# ============================================================

MODEL_FILE = "real_dna_cnn_model.keras"
DATA_FILE = "data/real_dna_dataset.csv"

OUTPUT_CSV = "real_motif_importance.csv"
OUTPUT_PNG = "real_motif_importance.png"


# ------------------------------------------------------------
# 1. Load trained CNN
# ------------------------------------------------------------

print("Loading trained CNN model...")

model = tf.keras.models.load_model(MODEL_FILE)

print("Model loaded successfully.")


# ------------------------------------------------------------
# 2. Load dataset
# ------------------------------------------------------------

df = pd.read_csv(DATA_FILE)

# Use only positive JASPAR CTCF sequences
positive_df = df[df["label"] == 1].copy()

print("\nTotal positive CTCF sequences:", len(positive_df))


# ------------------------------------------------------------
# 3. DNA one-hot encoding
# ------------------------------------------------------------

BASES = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}


def one_hot_encode(sequence):
    """
    Convert DNA sequence into 4-channel one-hot encoding.

    A = [1,0,0,0]
    C = [0,1,0,0]
    G = [0,0,1,0]
    T = [0,0,0,1]

    N = [0,0,0,0]
    """

    sequence = sequence.upper()

    encoded = np.zeros((len(sequence), 4), dtype=np.float32)

    for i, base in enumerate(sequence):
        if base in BASES:
            encoded[i, BASES[base]] = 1.0

    return encoded


# ------------------------------------------------------------
# 4. Calculate Gradient x Input saliency
# ------------------------------------------------------------

print("\nCalculating Gradient x Input saliency...")


all_importance = []


for index, row in positive_df.iterrows():

    sequence = str(row["sequence"])

    x = one_hot_encode(sequence)

    # Add batch dimension
    x_tensor = tf.convert_to_tensor(
        x[np.newaxis, :, :],
        dtype=tf.float32
    )

    with tf.GradientTape() as tape:

        tape.watch(x_tensor)

        prediction = model(x_tensor, training=False)

        # Probability of positive CTCF class
        score = prediction[:, 0]

    gradients = tape.gradient(score, x_tensor)

    # Gradient x Input
    saliency = gradients * x_tensor

    # Sum absolute importance across A,C,G,T
    position_importance = np.sum(
        np.abs(saliency.numpy()[0]),
        axis=1
    )

    all_importance.append(position_importance)


# ------------------------------------------------------------
# 5. Average importance across all CTCF sequences
# ------------------------------------------------------------

importance_matrix = np.array(all_importance)

mean_importance = np.mean(
    importance_matrix,
    axis=0
)

print("Importance matrix shape:", importance_matrix.shape)


# ------------------------------------------------------------
# 6. Save position-wise importance
# ------------------------------------------------------------

positions = np.arange(
    1,
    len(mean_importance) + 1
)

importance_df = pd.DataFrame({
    "position": positions,
    "importance": mean_importance
})

importance_df.to_csv(
    OUTPUT_CSV,
    index=False
)

print("\nSaved:", OUTPUT_CSV)


# ------------------------------------------------------------
# 7. Find top important positions
# ------------------------------------------------------------

top_n = min(20, len(mean_importance))

top_positions = (
    importance_df
    .sort_values("importance", ascending=False)
    .head(top_n)
)

print("\n======================================")
print("TOP IMPORTANT DNA POSITIONS")
print("======================================")

print(top_positions.to_string(index=False))


# ------------------------------------------------------------
# 8. Identify important region
# ------------------------------------------------------------

threshold = np.percentile(
    mean_importance,
    90
)

important_positions = positions[
    mean_importance >= threshold
]

print("\n90th percentile importance threshold:",
      round(float(threshold), 6))

print("\nImportant positions:")
print(important_positions.tolist())


# ------------------------------------------------------------
# 9. Plot motif importance
# ------------------------------------------------------------

plt.figure(figsize=(14, 5))

plt.plot(
    positions,
    mean_importance,
    linewidth=2
)

plt.xlabel("DNA sequence position")
plt.ylabel("Mean Gradient × Input importance")
plt.title("CNN Saliency Map for CTCF Motif Discovery")

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    OUTPUT_PNG,
    dpi=300
)

plt.close()

print("\nSaved:", OUTPUT_PNG)


# ------------------------------------------------------------
# 10. Print highest-importance region
# ------------------------------------------------------------

window_size = 10

if len(mean_importance) >= window_size:

    window_scores = np.convolve(
        mean_importance,
        np.ones(window_size),
        mode="valid"
    )

    best_start = int(
        np.argmax(window_scores)
    )

    best_end = best_start + window_size

    print("\n======================================")
    print("MOST IMPORTANT 10-BP REGION")
    print("======================================")

    print(
        "Positions:",
        best_start + 1,
        "to",
        best_end
    )

    print(
        "Mean importance:",
        round(
            float(window_scores[best_start] / window_size),
            6
        )
    )


print("\n======================================")
print("MOTIF INTERPRETABILITY COMPLETED")
print("======================================")
