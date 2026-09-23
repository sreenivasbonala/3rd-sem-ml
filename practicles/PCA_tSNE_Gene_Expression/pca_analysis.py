import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Load expression matrix
X = pd.read_csv("expression_matrix.csv", index_col=0)

# Transpose: samples × genes
X = X.T

print("Samples × Features:", X.shape)

# Standardize gene expression
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("Standardized matrix shape:", X_scaled.shape)

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

print("PCA shape:", X_pca.shape)
print("PC1 variance:", pca.explained_variance_ratio_[0])
print("PC2 variance:", pca.explained_variance_ratio_[1])
print("Total variance explained:", pca.explained_variance_ratio_.sum())

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.scatter(X_pca[:, 0], X_pca[:, 1], s=50)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
plt.title("PCA of GSE10072 Gene Expression")
plt.grid(True)
plt.tight_layout()

plt.savefig("PCA_plot.png", dpi=300)
plt.show()

print("Saved: PCA_plot.png")
