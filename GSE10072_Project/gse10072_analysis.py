# -*- coding: utf-8 -*-

# ==========================================
# GSE10072 GENE EXPRESSION ANALYSIS
# ==========================================

import GEOparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ==========================================
# STEP 1: LOAD GSE10072
# ==========================================

gse = GEOparse.get_GEO(
    geo="GSE10072",
    destdir="."
)

print("\nGSE10072 retrieved successfully!")

print("Number of samples:", len(gse.gsms))
print("Number of platforms:", len(gse.gpls))


# ==========================================
# STEP 2: CREATE EXPRESSION MATRIX
# ==========================================

expression = gse.pivot_samples("VALUE")

X = expression.copy()

print("\nExpression Matrix:")
print(expression.head())

print("\nOriginal expression shape:")
print(expression.shape)


# ==========================================
# STEP 3: CHECK MISSING VALUES
# ==========================================

missing = X.isnull().sum()

print("\nMissing values:")
print(missing)

print("\nTotal missing values:")
print(X.isnull().sum().sum())



# ==========================================
# STEP 4: REMOVE FEATURES WITH >10% MISSING
# ==========================================

X = X.dropna(
    axis=0,
    thresh=int(0.9 * X.shape[1])
)

print(
    "\nShape after removing features with >10% missing values:"
)

print(X.shape)


# ==========================================
# STEP 5: FILL REMAINING MISSING VALUES
# ==========================================

X = X.fillna(X.median())

print("\nTotal missing values after imputation:")

print(X.isnull().sum().sum())


# ==========================================
# STEP 6: TRANSPOSE EXPRESSION MATRIX
# ==========================================

expression_T = X.T

print("\n===================================")
print("TRANSPOSE")
print("===================================")

print("Original shape:")
print(X.shape)

print("\nTransposed shape:")
print(expression_T.shape)

print("\nTransposed expression matrix:")
print(expression_T.head())


# ==========================================
# STEP 7: CHECK EXPRESSION RANGE
# ==========================================

print("\n===================================")
print("EXPRESSION RANGE")
print("===================================")

print("\nMinimum expression value:")
print(X.min().min())

print("\nMaximum expression value:")
print(X.max().max())


# ==========================================
# STEP 8: KEEP ORIGINAL DATA
# ==========================================

X_original = expression_T.copy()

print("\nX_original shape:")
print(X_original.shape)


# ==========================================
# STEP 9: LOG2 TRANSFORMATION
# ==========================================

X_log = np.log2(X_original + 1)

print("\n===================================")
print("LOG2 TRANSFORMATION")
print("===================================")

print("\nOriginal expression values:")
print(X_original.iloc[:5, :5])

print("\nAfter Log2 transformation:")
print(X_log.iloc[:5, :5])

print("\nX_log shape:")
print(X_log.shape)


# ==========================================
# STEP 10: BEFORE LOG TRANSFORMATION PLOT
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(
    X_original.iloc[:, 0],
    bins=50
)

plt.xlabel("Expression")
plt.ylabel("Frequency")
plt.title("Expression Before Log Transformation")

plt.savefig(
    "Expression_Before_Log.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Before-log plot saved successfully!")


# ==========================================
# STEP 11: AFTER LOG TRANSFORMATION PLOT
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(
    X_log.iloc[:, 0],
    bins=50
)

plt.xlabel("log2(Expression + 1)")
plt.ylabel("Frequency")
plt.title("Expression After Log Transformation")

plt.savefig(
    "Expression_After_Log.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("After-log plot saved successfully!")


# ==========================================
# STEP 12: FEATURE VARIANCE
# ==========================================

feature_variance = X_log.var(axis=0)

print("\n===================================")
print("FEATURE VARIANCE")
print("===================================")

print("\nFeature variance:")
print(feature_variance.head())

print("\nFeature variance statistics:")
print(feature_variance.describe())


# ==========================================
# STEP 13: SELECT TOP 5000 FEATURES
# ==========================================

top_features = feature_variance.sort_values(
    ascending=False
).head(5000).index

print("\n===================================")
print("TOP 5000 FEATURES")
print("===================================")

print("\nNumber of top features:")
print(len(top_features))

print("\nTop 20 features:")
print(top_features[:20])


# ==========================================
# STEP 14: CREATE HIGH VARIANCE FEATURE MATRIX
# ==========================================

X_hvg = X_log[top_features]

print("\n===================================")
print("HIGH VARIANCE FEATURES")
print("===================================")

print("\nOriginal X_log shape:")
print(X_log.shape)

print("\nX_hvg shape:")
print(X_hvg.shape)

print("\nX_hvg first 5 rows:")
print(X_hvg.iloc[:5, :5])


# ==========================================
# STEP 15: IMPORT SKLEARN MODULES
# ==========================================

from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA

print("\n===================================")
print("SKLEARN MODULES")
print("===================================")

print("StandardScaler, VarianceThreshold and PCA imported successfully!")


# ==========================================
# STEP 16: STANDARDIZATION
# ==========================================

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X_hvg)

print("\n===================================")
print("STANDARDIZATION")
print("===================================")

print("X_hvg shape:")
print(X_hvg.shape)

print("X_scaled shape:")
print(X_scaled.shape)


# ==========================================
# CONVERT X_scaled TO DATAFRAME
# ==========================================

X_scaled = pd.DataFrame(
    X_scaled,
    index=X_hvg.index,
    columns=X_hvg.columns
)

print("\nX_scaled DataFrame:")
print(X_scaled.head())


# ==========================================

# CHECK MEAN AND STANDARD DEVIATION

# ==========================================

print("\n===================================")

print("MEAN AND STANDARD DEVIATION")

print("===================================")

print("\nMean of each feature:")

print(X_scaled.mean(axis=0).head())

print("\nStandard deviation of each feature:")

print(X_scaled.std(axis=0, ddof=0).head())

print("\nOverall mean:")

print(X_scaled.mean().mean())

print("\nOverall standard deviation:")

print(X_scaled.values.std())

# ==========================================
# STEP 17: VARIANCE THRESHOLD
# ==========================================

selector = VarianceThreshold(
    threshold=0.0
)

X_selected = selector.fit_transform(X_scaled)

print("\n===================================")
print("VARIANCE THRESHOLD")
print("===================================")

print("\nBefore VarianceThreshold:")
print(X_scaled.shape)

print("\nAfter VarianceThreshold:")
print(X_selected.shape)


# ==========================================
# STEP 18: PCA
# ==========================================

pca = PCA(
    n_components=2
)

X_pca = pca.fit_transform(
    X_selected
)

print("\n===================================")
print("PCA")
print("===================================")

print("\nPCA shape:")
print(X_pca.shape)

print("\nExplained variance ratio:")
print(pca.explained_variance_ratio_)

print("\nTotal explained variance:")
print(
    pca.explained_variance_ratio_.sum()
)
# ==========================================
# STEP 19: PCA PLOT
# ==========================================

print("\n===================================")
print("PCA PLOT")
print("===================================")

print("PCA plot started")
print("X_pca shape:", X_pca.shape)

plt.figure(figsize=(8, 5))

plt.scatter(
    X_pca[:, 0],
    X_pca[:, 1]
)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA of GSE10072 High-Variance Features")

plt.savefig(
    "GSE10072_PCA.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("PCA plot saved successfully!")
print("File: GSE10072_PCA.png")


# ==========================================
# STEP 20: FINAL SUMMARY
# ==========================================

print("\n===================================")
print("FINAL ANALYSIS SUMMARY")
print("===================================")

print("Original expression matrix:")
print(X.shape)

print("Log transformed matrix:")
print(X_log.shape)

print("Top 5000 feature matrix:")
print(X_hvg.shape)

print("Standardized matrix:")
print(X_scaled.shape)

print("Variance threshold matrix: ")
print(X_selected.shape)

print("PCA matrix:")
print(X_pca.shape)

print("\n===================================")
print("Analysis completed successfully!")
print("===================================")