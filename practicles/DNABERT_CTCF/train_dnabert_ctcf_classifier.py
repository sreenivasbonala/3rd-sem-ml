import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    classification_report
)

# --------------------------------------------------
# 1. Load DNABERT-2 embeddings
# --------------------------------------------------

positive_embeddings = np.load(
    "ctcf_dnabert2_embeddings.npy"
)

negative_embeddings = np.load(
    "ctcf_negative_dnabert2_embeddings.npy"
)

print("Positive embeddings:", positive_embeddings.shape)
print("Negative embeddings:", negative_embeddings.shape)

# --------------------------------------------------
# 2. Create labels
# --------------------------------------------------

positive_labels = np.ones(len(positive_embeddings))
negative_labels = np.zeros(len(negative_embeddings))

# --------------------------------------------------
# 3. Combine data
# --------------------------------------------------

X = np.vstack([
    positive_embeddings,
    negative_embeddings
])

y = np.concatenate([
    positive_labels,
    negative_labels
])

print("Combined feature matrix:", X.shape)
print("Combined labels:", y.shape)

# --------------------------------------------------
# 4. Train-test split
# --------------------------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("Training samples:", X_train.shape[0])
print("Testing samples:", X_test.shape[0])

# --------------------------------------------------
# 5. Feature scaling
# --------------------------------------------------

scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# --------------------------------------------------
# 6. Train Logistic Regression
# --------------------------------------------------

model = LogisticRegression(
    max_iter=2000,
    random_state=42
)

print("Training Logistic Regression...")

model.fit(X_train_scaled, y_train)

print("Training completed!")

# --------------------------------------------------
# 7. Predictions
# --------------------------------------------------

y_pred = model.predict(X_test_scaled)

y_prob = model.predict_proba(X_test_scaled)[:, 1]

# --------------------------------------------------
# 8. Evaluation
# --------------------------------------------------

accuracy = accuracy_score(y_test, y_pred)

precision = precision_score(y_test, y_pred)

recall = recall_score(y_test, y_pred)

f1 = f1_score(y_test, y_pred)

auc = roc_auc_score(y_test, y_prob)

cm = confusion_matrix(y_test, y_pred)

print("\n========== DNABERT-2 CTCF CLASSIFIER ==========")

print("Accuracy :", round(accuracy, 4))
print("Precision:", round(precision, 4))
print("Recall   :", round(recall, 4))
print("F1-score :", round(f1, 4))
print("ROC-AUC  :", round(auc, 4))

print("\nConfusion Matrix:")
print(cm)

print("\nClassification Report:")
print(
    classification_report(
        y_test,
        y_pred,
        target_names=["Negative", "Positive"]
    )
)

# --------------------------------------------------
# 9. Save predictions
# --------------------------------------------------

results = pd.DataFrame({
    "Actual": y_test.astype(int),
    "Predicted": y_pred.astype(int),
    "Probability_Positive": y_prob
})

results.to_csv(
    "dnabert_ctcf_test_predictions.csv",
    index=False
)

# --------------------------------------------------
# 10. Save combined embeddings
# --------------------------------------------------

np.save(
    "ctcf_all_dnabert2_embeddings.npy",
    X
)

# Save labels
pd.DataFrame({
    "label": y.astype(int)
}).to_csv(
    "ctcf_labels.csv",
    index=False
)

print("\nSaved:")
print("  ctcf_all_dnabert2_embeddings.npy")
print("  ctcf_labels.csv")
print("  dnabert_ctcf_test_predictions.csv")

print("\nDNABERT-2 CTCF classification completed!")
