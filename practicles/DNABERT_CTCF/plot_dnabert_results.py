import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    confusion_matrix,
    ConfusionMatrixDisplay,
    roc_curve,
    auc
)

# Load predictions
df = pd.read_csv("dnabert_ctcf_test_predictions.csv")

print("Columns:", df.columns.tolist())
print("Number of test samples:", len(df))

# Correct column names
y_true = df["Actual"]
y_pred = df["Predicted"]
y_prob = df["Probability_Positive"]

# --------------------------------------------------
# 1. Confusion Matrix
# --------------------------------------------------

cm = confusion_matrix(y_true, y_pred)

print("\nConfusion Matrix:")
print(cm)

disp = ConfusionMatrixDisplay(
    confusion_matrix=cm,
    display_labels=["Negative", "Positive"]
)

disp.plot()
plt.title("DNABERT-2 CTCF Classification - Confusion Matrix")
plt.tight_layout()

plt.savefig(
    "dnabert_ctcf_confusion_matrix.png",
    dpi=300
)

plt.close()

# --------------------------------------------------
# 2. ROC Curve
# --------------------------------------------------

fpr, tpr, thresholds = roc_curve(y_true, y_prob)

roc_auc = auc(fpr, tpr)

print("\nROC-AUC:", round(roc_auc, 4))

plt.figure(figsize=(7, 6))

plt.plot(
    fpr,
    tpr,
    label=f"DNABERT-2 (AUC = {roc_auc:.4f})"
)

plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--",
    label="Random classifier"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("DNABERT-2 CTCF Classification - ROC Curve")

plt.legend()
plt.grid(True)

plt.tight_layout()

plt.savefig(
    "dnabert_ctcf_roc_curve.png",
    dpi=300
)

plt.close()

# --------------------------------------------------
# 3. Summary
# --------------------------------------------------

print("\nSaved:")
print("  dnabert_ctcf_confusion_matrix.png")
print("  dnabert_ctcf_roc_curve.png")

print("\nDNABERT-2 result plotting completed!")
