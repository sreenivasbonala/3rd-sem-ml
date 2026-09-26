import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

# ==============================
# Load dataset and model
# ==============================

df = pd.read_csv("dna_motif_dataset.csv")
model = tf.keras.models.load_model("cnn_motif_model.keras")

mapping = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1]
}


def one_hot_encode(sequence):
    return np.array(
        [mapping[base] for base in sequence],
        dtype=np.float32
    )


# ==============================
# Select a positive sequence
# ==============================

positive_df = df[df["label"] == 1]

sequence = positive_df.iloc[0]["sequence"]

X = one_hot_encode(sequence)
X = np.expand_dims(X, axis=0)


# ==============================
# Prediction
# ==============================

prediction = model.predict(X, verbose=0)[0][0]

print("====================================")
print("CNN MOTIF INTERPRETABILITY")
print("====================================")

print("Sequence:")
print(sequence)

print("\nKnown motif:")
print("CACGTG")

print("\nCNN prediction probability:")
print(round(float(prediction), 4))


# ==============================
# Gradient-based saliency
# ==============================

input_tensor = tf.convert_to_tensor(X)

with tf.GradientTape() as tape:

    tape.watch(input_tensor)

    prediction_value = model(
        input_tensor,
        training=False
    )[0, 0]

gradients = tape.gradient(
    prediction_value,
    input_tensor
).numpy()[0]


# ==============================
# Calculate importance
# ==============================

importance = np.abs(
    gradients * X[0]
).sum(axis=1)


# ==============================
# Find most important 6-mer
# ==============================

window_size = 6

window_scores = []

for i in range(len(sequence) - window_size + 1):

    score = importance[
        i:i + window_size
    ].sum()

    window_scores.append(score)

best_start = int(
    np.argmax(window_scores)
)

best_end = best_start + window_size

important_motif = sequence[
    best_start:best_end
]

print("\nMost important region:")
print(
    f"Positions {best_start + 1}-{best_end}"
)

print("Detected 6-mer:", important_motif)

print(
    "\nExpected motif: CACGTG"
)


# ==============================
# Save importance table
# ==============================

importance_df = pd.DataFrame({
    "Position": np.arange(1, len(sequence) + 1),
    "Base": list(sequence),
    "Importance": importance
})

importance_df.to_csv(
    "motif_importance.csv",
    index=False
)


# ==============================
# Plot motif importance
# ==============================

plt.figure(figsize=(12, 5))

plt.bar(
    np.arange(1, len(sequence) + 1),
    importance
)

plt.xlabel("DNA Position")
plt.ylabel("Importance")
plt.title("CNN Gradient-Based DNA Motif Importance")

plt.xticks(
    np.arange(1, len(sequence) + 1, 5)
)

plt.tight_layout()

plt.savefig(
    "motif_importance.png",
    dpi=300
)

plt.close()


print("\nSaved:")
print("motif_importance.csv")
print("motif_importance.png")

print("\nInterpretability analysis completed.")
