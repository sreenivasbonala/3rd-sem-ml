import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Load expression data
X = pd.read_csv("expression_matrix.csv", index_col=0)
X = X.T

# Load metadata
metadata = pd.read_csv("sample_metadata.csv")

# Make sure sample order matches
metadata = metadata.set_index("Sample")
metadata = metadata.loc[X.index]

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# =========================
# PCA
# =========================

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

plt.figure(figsize=(8, 6))

for group in ["Lung Tumor", "Normal Lung"]:

    mask = metadata["Group"] == group

    plt.scatter(
        X_pca[mask, 0],
        X_pca[mask, 1],
        label=group,
        s=50
    )

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
plt.title("PCA of GSE10072: Lung Tumor vs Normal Lung")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("PCA_Tumor_Normal.png", dpi=300)

# =========================
# t-SNE
# =========================

# PCA to 30 dimensions before t-SNE
pca30 = PCA(n_components=30)
X_pca30 = pca30.fit_transform(X_scaled)

tsne = TSNE(
    n_components=2,
    perplexity=30,
    random_state=42,
    init="pca",
    learning_rate="auto",
    max_iter=1000
)

X_tsne = tsne.fit_transform(X_pca30)

plt.figure(figsize=(8, 6))

for group in ["Lung Tumor", "Normal Lung"]:

    mask = metadata["Group"] == group

    plt.scatter(
        X_tsne[mask, 0],
        X_tsne[mask, 1],
        label=group,
        s=50
    )

plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.title("t-SNE of GSE10072: Lung Tumor vs Normal Lung")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("tSNE_Tumor_Normal.png", dpi=300)

print("Saved: PCA_Tumor_Normal.png")
print("Saved: tSNE_Tumor_Normal.png")

