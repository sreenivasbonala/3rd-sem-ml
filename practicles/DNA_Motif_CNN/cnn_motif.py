import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve

# ==============================
# 1. Load dataset
# ==============================

df = pd.read_csv("dna_motif_dataset.csv")

sequences = df["sequence"].values
labels = df["label"].values

print("Dataset shape:", df.shape)
print("Positive samples:", sum(labels == 1))
print("Negative samples:", sum(labels == 0))


# ==============================
# 2. One-hot encoding
# ==============================

mapping = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1]
}


def one_hot_encode(sequence):
    return np.array([mapping[base] for base in sequence], dtype=np.float32)


X = np.array([one_hot_encode(seq) for seq in sequences])
y = np.array(labels)

print("Encoded data shape:", X.shape)


# ==============================
# 3. Train-test split
# ==============================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("Training samples:", len(X_train))
print("Testing samples:", len(X_test))


# ==============================
# 4. Build CNN
# ==============================

model = tf.keras.Sequential([
    tf.keras.layers.Input(shape=(50, 4)),

    tf.keras.layers.Conv1D(
        filters=32,
        kernel_size=6,
        activation="relu"
    ),

    tf.keras.layers.GlobalMaxPooling1D(),

    tf.keras.layers.Dense(
        32,
        activation="relu"
    ),

    tf.keras.layers.Dropout(0.3),

    tf.keras.layers.Dense(
        1,
        activation="sigmoid"
    )
])

model.compile(
    optimizer="adam",
    loss="binary_crossentropy",
    metrics=[
        tf.keras.metrics.AUC(name="auc")
    ]
)

print("\nCNN Model:")
model.summary()


# ==============================
# 5. Train CNN
# ==============================

early_stop = tf.keras.callbacks.EarlyStopping(
    monitor="val_auc",
    mode="max",
    patience=5,
    restore_best_weights=True
)

history = model.fit(
    X_train,
    y_train,
    validation_split=0.20,
    epochs=20,
    batch_size=64,
    callbacks=[early_stop],
    verbose=1
)


# ==============================
# 6. Test prediction
# ==============================

y_prob = model.predict(X_test).ravel()

test_auc = roc_auc_score(y_test, y_prob)

print("\n====================================")
print("CNN MOTIF DISCOVERY RESULTS")
print("====================================")
print("Test AUROC:", round(test_auc, 4))


# ==============================
# 7. Save predictions
# ==============================

results = pd.DataFrame({
    "sequence": [
        "".join(
            "ACGT"[np.argmax(base)]
            for base in seq
        )
        for seq in X_test
    ],
    "actual_label": y_test,
    "predicted_probability": y_prob
})

results.to_csv(
    "test_predictions.csv",
    index=False
)


# ==============================
# 8. Training curve
# ==============================

plt.figure(figsize=(8, 5))

plt.plot(
    history.history["auc"],
    label="Training AUC"
)

plt.plot(
    history.history["val_auc"],
    label="Validation AUC"
)

plt.xlabel("Epoch")
plt.ylabel("AUROC")
plt.title("CNN Training and Validation AUROC")
plt.legend()
plt.tight_layout()

plt.savefig(
    "training_curve.png",
    dpi=300
)

plt.close()


# ==============================
# 9. ROC curve
# ==============================

fpr, tpr, thresholds = roc_curve(
    y_test,
    y_prob
)

plt.figure(figsize=(7, 6))

plt.plot(
    fpr,
    tpr,
    label=f"AUROC = {test_auc:.4f}"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("CNN DNA Motif Detection ROC Curve")
plt.legend()
plt.tight_layout()

plt.savefig(
    "roc_curve.png",
    dpi=300
)

plt.close()


# ==============================
# 10. Save model
# ==============================

model.save("cnn_motif_model.keras")

print("\nSaved files:")
print("cnn_motif_model.keras")
print("training_curve.png")
print("roc_curve.png")
print("test_predictions.csv")

print("\nAnalysis completed successfully.")
