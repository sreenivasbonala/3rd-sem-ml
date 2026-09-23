import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Load expression matrix
X = pd.read_csv("expression_matrix.csv", index_col=0)

# Samples × genes
X = X.T

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA to reduce dimensions before t-SNE
pca = PCA(n_components=30)
X_pca = pca.fit_transform(X_scaled)

print("PCA-reduced shape:", X_pca.shape)

# t-SNE
tsne = TSNE(
    n_components=2,
    perplexity=30,
    random_state=42,
    init="pca",
    learning_rate="auto",
    max_iter=1000
)

X_tsne = tsne.fit_transform(X_pca)

print("t-SNE shape:", X_tsne.shape)

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(X_tsne[:, 0], X_tsne[:, 1], s=50)

plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.title("t-SNE of GSE10072 Gene Expression")
plt.grid(True)
plt.tight_layout()

plt.savefig("tSNE_plot.png", dpi=300)
plt.show()

print("Saved: tSNE_plot.png")
