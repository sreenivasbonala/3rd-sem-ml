import os
import random
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    classification_report,
    confusion_matrix,
    roc_curve,
    precision_recall_curve
)

# ============================================================
# 1. Reproducibility
# ============================================================

SEED = 42

random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)

print("TensorFlow version:", tf.__version__)

# ============================================================
# 2. Load real biological dataset
# ============================================================

DATA_FILE = "data/real_dna_dataset.csv"

df = pd.read_csv(DATA_FILE)

print("\nDataset loaded")
print("Total sequences:", len(df))
print("Sequence length:", df["sequence"].str.len().unique())
print("\nClass distribution:")
print(df["label"].value_counts())

# ============================================================
# 3. DNA one-hot encoding
# ============================================================

BASES = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
    "N": [0, 0, 0, 0]
}


def one_hot_encode(sequence):
    return np.array(
        [BASES.get(base, [0, 0, 0, 0]) for base in sequence],
        dtype=np.float32
    )


X = np.array([
    one_hot_encode(seq)
    for seq in df["sequence"]
])

y = df["label"].values.astype(np.float32)

print("\nEncoded data shape:", X.shape)
print("Labels shape:", y.shape)

# ============================================================
# 4. Train / validation / test split
# ============================================================

indices = np.arange(len(df))

train_idx, temp_idx = train_test_split(
    indices,
    test_size=0.30,
    stratify=y,
    random_state=SEED
)

val_idx, test_idx = train_test_split(
    temp_idx,
    test_size=0.50,
    stratify=y[temp_idx],
    random_state=SEED
)

X_train = X[train_idx]
X_val = X[val_idx]
X_test = X[test_idx]

y_train = y[train_idx]
y_val = y[val_idx]
y_test = y[test_idx]

print("\nData split:")
print("Training:", len(X_train))
print("Validation:", len(X_val))
print("Testing:", len(X_test))

print("\nTraining class distribution:")
print(pd.Series(y_train).value_counts())

print("\nValidation class distribution:")
print(pd.Series(y_val).value_counts())

print("\nTesting class distribution:")
print(pd.Series(y_test).value_counts())

# ============================================================
# 5. Build CNN model
# ============================================================

model = tf.keras.Sequential([
    tf.keras.layers.Input(shape=(201, 4)),

    tf.keras.layers.Conv1D(
        filters=64,
        kernel_size=8,
        activation="relu"
    ),

    tf.keras.layers.GlobalMaxPooling1D(),

    tf.keras.layers.Dense(
        32,
        activation="relu"
    ),

    tf.keras.layers.Dropout(0.30),

    tf.keras.layers.Dense(
        1,
        activation="sigmoid"
    )
])

model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
    loss="binary_crossentropy",
    metrics=[
        "accuracy",
        tf.keras.metrics.AUC(name="auc"),
        tf.keras.metrics.AUC(
            name="auprc",
            curve="PR"
        )
    ]
)

print("\nCNN model:")
model.summary()

# ============================================================
# 6. Early stopping
# ============================================================

early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor="val_auc",
    mode="max",
    patience=8,
    restore_best_weights=True
)

# ============================================================
# 7. Train CNN
# ============================================================

print("\nStarting CNN training...")

history = model.fit(
    X_train,
    y_train,
    validation_data=(X_val, y_val),
    epochs=50,
    batch_size=32,
    callbacks=[early_stopping],
    verbose=1
)

# ============================================================
# 8. Save trained model
# ============================================================

model.save("real_dna_cnn_model.keras")

print("\nModel saved as: real_dna_cnn_model.keras")

# ============================================================
# 9. Test predictions
# ============================================================

y_prob = model.predict(
    X_test,
    verbose=0
).ravel()

y_pred = (y_prob >= 0.5).astype(int)

# ============================================================
# 10. Evaluation metrics
# ============================================================

accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred, zero_division=0)
recall = recall_score(y_test, y_pred, zero_division=0)
f1 = f1_score(y_test, y_pred, zero_division=0)

auroc = roc_auc_score(y_test, y_prob)
auprc = average_precision_score(y_test, y_prob)

print("\n======================================")
print("REAL CTCF CNN TEST RESULTS")
print("======================================")

print(f"Accuracy : {accuracy:.4f}")
print(f"Precision: {precision:.4f}")
print(f"Recall   : {recall:.4f}")
print(f"F1-score : {f1:.4f}")
print(f"AUROC    : {auroc:.4f}")
print(f"AUPRC    : {auprc:.4f}")

print("\nClassification report:")
print(
    classification_report(
        y_test,
        y_pred,
        target_names=["Shuffled Control", "JASPAR CTCF"],
        zero_division=0
    )
)

print("\nConfusion matrix:")
print(confusion_matrix(y_test, y_pred))

# ============================================================
# 11. Save test predictions
# ============================================================

results = df.iloc[test_idx].copy()

results["predicted_probability"] = y_prob
results["predicted_label"] = y_pred

results.to_csv(
    "real_test_predictions.csv",
    index=False
)

print("\nSaved: real_test_predictions.csv")

# ============================================================
# 12. Training curves
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(
    history.history["loss"],
    label="Training Loss"
)

plt.plot(
    history.history["val_loss"],
    label="Validation Loss"
)

plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Real CTCF CNN Training and Validation Loss")
plt.legend()
plt.tight_layout()

plt.savefig(
    "real_training_curve.png",
    dpi=300
)

plt.close()

print("Saved: real_training_curve.png")

# ============================================================
# 13. ROC curve
# ============================================================

fpr, tpr, _ = roc_curve(y_test, y_prob)

plt.figure(figsize=(6, 6))

plt.plot(
    fpr,
    tpr,
    label=f"AUROC = {auroc:.4f}"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Real CTCF CNN ROC Curve")
plt.legend()
plt.tight_layout()

plt.savefig(
    "real_roc_curve.png",
    dpi=300
)

plt.close()

print("Saved: real_roc_curve.png")

# ============================================================
# 14. Precision-Recall curve
# ============================================================

precision_curve, recall_curve, _ = precision_recall_curve(
    y_test,
    y_prob
)

plt.figure(figsize=(6, 6))

plt.plot(
    recall_curve,
    precision_curve,
    label=f"AUPRC = {auprc:.4f}"
)

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Real CTCF CNN Precision-Recall Curve")
plt.legend()
plt.tight_layout()

plt.savefig(
    "real_pr_curve.png",
    dpi=300
)

plt.close()

print("Saved: real_pr_curve.png")

# ============================================================
# 15. Final summary
# ============================================================

print("\n======================================")
print("FILES GENERATED")
print("======================================")

print("real_dna_cnn_model.keras")
print("real_test_predictions.csv")
print("real_training_curve.png")
print("real_roc_curve.png")
print("real_pr_curve.png")

print("\nExperiment completed successfully.")


