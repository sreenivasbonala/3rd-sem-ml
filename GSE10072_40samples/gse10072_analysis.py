import GEOparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif


# ==========================================
# STEP 1: LOAD GSE10072
# ==========================================

gse = GEOparse.get_GEO(
    geo="GSE10072",
    destdir="."
)

print("\nGSE10072 loaded successfully")
print("Total samples:", len(gse.gsms))


# ==========================================
# STEP 2: IDENTIFY TUMOR AND NORMAL SAMPLES
# ==========================================

tumor_samples = []
normal_samples = []

for gsm_name, gsm in gse.gsms.items():

    title = gsm.metadata.get("title", [""])[0]

    if "Lung Tumor" in title:
        tumor_samples.append(gsm_name)

    elif "Normal Lung" in title:
        normal_samples.append(gsm_name)


print("\nOriginal group counts:")
print("Lung Tumor:", len(tumor_samples))
print("Normal Lung:", len(normal_samples))


# ==========================================
# SELECT 20 TUMOR + 20 NORMAL
# ==========================================

tumor_20 = tumor_samples[:20]
normal_20 = normal_samples[:20]

selected_samples = tumor_20 + normal_20


print("\n================================")
print("SELECTED 40 SAMPLES")
print("================================")

print("\nLung Tumor samples:")

for sample in tumor_20:
    print(sample)

print("\nNormal Lung samples:")

for sample in normal_20:
    print(sample)

print("\nTotal selected samples:",
      len(selected_samples))


# ==========================================
# STEP 3: CREATE EXPRESSION MATRIX
# ==========================================

expression_data = []

for gsm_name in selected_samples:

    gsm = gse.gsms[gsm_name]

    table = gsm.table

    sample_data = table[["ID_REF", "VALUE"]].copy()

    sample_data = sample_data.rename(
        columns={"VALUE": gsm_name}
    )

    sample_data = sample_data.set_index("ID_REF")

    expression_data.append(sample_data)


expression_matrix = pd.concat(
    expression_data,
    axis=1
)


print("\n================================")
print("ORIGINAL EXPRESSION MATRIX")
print("================================")

print("Rows:", expression_matrix.shape[0])
print("Columns:", expression_matrix.shape[1])

print("Matrix shape:",
      expression_matrix.shape)

print("\nFirst 5 rows:")
print(expression_matrix.head())


# ==========================================
# TRANSPOSE
# ==========================================

X = expression_matrix.T


print("\n================================")
print("TRANSPOSED EXPRESSION MATRIX")
print("================================")

print("Rows (samples):", X.shape[0])
print("Columns (features):", X.shape[1])

print("Matrix shape:",
      X.shape)


# ==========================================
# SAVE EXPRESSION MATRIX
# ==========================================

X.to_csv(
    "GSE10072_40samples_expression_matrix.csv"
)

print("\nExpression matrix saved successfully!")


# ==========================================
# STEP 4: MISSING VALUES
# ==========================================

print("\n================================")
print("MISSING VALUE CHECK")
print("================================")

missing_values = X.isnull().sum().sum()

print("Total missing values:",
      missing_values)


if missing_values > 0:

    print("Missing values found.")
    print("Filling with feature median.")

    X = X.fillna(X.median())

else:

    print("No missing values found.")
    print("No missing value treatment required.")


# ==========================================
# DUPLICATE FEATURES
# ==========================================

print("\n================================")
print("DUPLICATE FEATURE CHECK")
print("================================")

duplicate_features = X.columns.duplicated().sum()

print("Number of duplicate features:",
      duplicate_features)


if duplicate_features > 0:

    X = X.loc[
        :,
        ~X.columns.duplicated()
    ]

    print("Duplicate features removed.")

else:

    print("No duplicate features found.")


# ==========================================
# CLEAN MATRIX
# ==========================================

print("\n================================")
print("MATRIX AFTER CLEANING")
print("================================")

print("Rows (samples):",
      X.shape[0])

print("Columns (features):",
      X.shape[1])


X.to_csv(
    "GSE10072_40samples_clean_matrix.csv"
)

print("\nClean expression matrix saved successfully!")


# ==========================================
# STEP 5: DATA DISTRIBUTION
# ==========================================

print("\n================================")
print("DATA DISTRIBUTION")
print("================================")

print("Minimum expression value:",
      X.min().min())

print("Maximum expression value:",
      X.max().max())

print("Mean expression value:",
      X.mean().mean())

print("Median expression value:",
      X.median().median())


# ==========================================
# HISTOGRAM
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(
    X.values.flatten(),
    bins=50
)

plt.xlabel("Expression Value")
plt.ylabel("Frequency")

plt.title(
    "GSE10072 Expression Distribution"
)

plt.tight_layout()

plt.savefig(
    "expression_distribution_before_normalization.png"
)

plt.close()

print("\nDistribution plot saved successfully!")


# ==========================================
# STANDARDIZATION
# ==========================================

print("\n================================")
print("NORMALIZATION / STANDARDIZATION")
print("================================")

scaler = StandardScaler()

X_scaled_array = scaler.fit_transform(X)


X_scaled = pd.DataFrame(
    X_scaled_array,
    index=X.index,
    columns=X.columns
)


print("\nAfter standardization:")

print(
    "Mean:",
    X_scaled.mean().mean()
)

print(
    "Standard deviation:",
    X_scaled.std().mean()
)


# ==========================================
# SAVE NORMALIZED MATRIX
# ==========================================

X_scaled.to_csv(
    "GSE10072_40samples_normalized.csv"
)

print("\nNormalized matrix saved successfully!")


# ==========================================
# STEP 6: FEATURE SELECTION
# ==========================================

print("\n================================")
print("FEATURE SELECTION")
print("================================")


# ==========================================
# CREATE CLASS LABELS
# ==========================================

# 0 = Normal Lung
# 1 = Lung Tumor

y = []

for sample in selected_samples:

    if sample in tumor_samples:

        y.append(1)

    elif sample in normal_samples:

        y.append(0)


y = np.array(y)


print("\nClass labels created")

print("Number of Normal Lung samples:",
      np.sum(y == 0))

print("Number of Lung Tumor samples:",
      np.sum(y == 1))


# ==========================================
# METHOD 1: VARIANCE SELECTION
# ==========================================

print("\n--------------------------------")
print("METHOD 1: VARIANCE SELECTION")
print("--------------------------------")


# Calculate variance
feature_variance = X_scaled.var(
    axis=0
)


# Select top 5000 features
top_5000_variance = (
    feature_variance
    .sort_values(ascending=False)
    .head(5000)
)


# Create selected matrix
variance_matrix = X_scaled[
    top_5000_variance.index
]


print(
    "Original number of features:",
    X_scaled.shape[1]
)

print(
    "Selected features:",
    variance_matrix.shape[1]
)

print(
    "Variance-selected matrix shape:",
    variance_matrix.shape
)


# Save
variance_matrix.to_csv(
    "GSE10072_40samples_variance_top5000.csv"
)

print(
    "Variance-selected matrix saved!"
)


# ==========================================
# METHOD 2: ANOVA F-TEST
# ==========================================

print("\n--------------------------------")
print("METHOD 2: ANOVA F-TEST")
print("--------------------------------")


selector = SelectKBest(
    score_func=f_classif,
    k=5000
)


X_anova_array = selector.fit_transform(
    X_scaled,
    y
)


# Get selected feature names
anova_features = X_scaled.columns[
    selector.get_support()
]


# Convert to DataFrame
anova_matrix = pd.DataFrame(
    X_anova_array,
    index=X_scaled.index,
    columns=anova_features
)


print(
    "Original number of features:",
    X_scaled.shape[1]
)

print(
    "Selected features:",
    anova_matrix.shape[1]
)

print(
    "ANOVA-selected matrix shape:",
    anova_matrix.shape
)


# Save
anova_matrix.to_csv(
    "GSE10072_40samples_ANOVA_top5000.csv"
)

print(
    "ANOVA-selected matrix saved!"
)


# ==========================================
# FEATURE SELECTION SUMMARY
# ==========================================

print("\n================================")
print("FEATURE SELECTION SUMMARY")
print("================================")

print(
    "Original features:",
    X_scaled.shape[1]
)

print(
    "Variance-selected features:",
    variance_matrix.shape[1]
)

print(
    "ANOVA-selected features:",
    anova_matrix.shape[1]
)


print("\nFirst 10 variance-selected features:")

print(
    list(variance_matrix.columns[:10])
)


print("\nFirst 10 ANOVA-selected features:")

print(
    list(anova_matrix.columns[:10])
)


print("\n================================")
print("STEP 6 COMPLETED")
print("================================")
# =================================================
# STEP 7: PCA
# =================================================

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

print("\n================================")
print("STEP 7: PCA ANALYSIS")
print("================================")


# =================================================
# PCA ON VARIANCE-SELECTED FEATURES
# =================================================

print("\n--------------------------------")
print("PCA: VARIANCE-SELECTED FEATURES")
print("--------------------------------")

pca_variance = PCA(n_components=2)

variance_pca = pca_variance.fit_transform(
    variance_matrix
)

print("PCA output shape:",
      variance_pca.shape)

print("Explained variance ratio:")
print(pca_variance.explained_variance_ratio_)

print("Total explained variance:",
      pca_variance.explained_variance_ratio_.sum())


# =================================================
# PCA PLOT - VARIANCE
# =================================================

plt.figure(figsize=(8, 6))

plt.scatter(
    variance_pca[y == 0, 0],
    variance_pca[y == 0, 1],
    label="Normal Lung"
)

plt.scatter(
    variance_pca[y == 1, 0],
    variance_pca[y == 1, 1],
    label="Lung Tumor"
)

plt.xlabel(
    "PC1"
)

plt.ylabel(
    "PC2"
)

plt.title(
    "PCA - Variance Selected Features"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "PCA_variance_selected.png",
    dpi=300
)

plt.close()

print("Variance PCA plot saved!")


# =================================================
# PCA ON ANOVA-SELECTED FEATURES
# =================================================

print("\n--------------------------------")
print("PCA: ANOVA-SELECTED FEATURES")
print("--------------------------------")

pca_anova = PCA(n_components=2)

anova_pca = pca_anova.fit_transform(
    anova_matrix
)

print("PCA output shape:",
      anova_pca.shape)

print("Explained variance ratio:")
print(pca_anova.explained_variance_ratio_)

print("Total explained variance:",
      pca_anova.explained_variance_ratio_.sum())


# =================================================
# PCA PLOT - ANOVA
# =================================================

plt.figure(figsize=(8, 6))

plt.scatter(
    anova_pca[y == 0, 0],
    anova_pca[y == 0, 1],
    label="Normal Lung"
)

plt.scatter(
    anova_pca[y == 1, 0],
    anova_pca[y == 1, 1],
    label="Lung Tumor"
)

plt.xlabel(
    "PC1"
)

plt.ylabel(
    "PC2"
)

plt.title(
    "PCA - ANOVA Selected Features"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    "PCA_ANOVA_selected.png",
    dpi=300
)

plt.close()

print("ANOVA PCA plot saved!")


# =================================================
# PCA SUMMARY
# =================================================

print("\n================================")
print("PCA SUMMARY")
print("================================")

print(
    "Variance PCA explained variance:",
    pca_variance.explained_variance_ratio_.sum()
)

print(
    "ANOVA PCA explained variance:",
    pca_anova.explained_variance_ratio_.sum()
)

print("\n================================")
print("STEP 7 COMPLETED")
print("================================")
# =================================================
# STEP 9: COMPARE FEATURE SELECTION METHODS
# =================================================

print("\n================================")
print("STEP 9: FEATURE SELECTION COMPARISON")
print("================================")


# =================================================
# 1. NUMBER OF FEATURES
# =================================================

print("\nNumber of selected features:")

print("Variance method:",
      variance_matrix.shape[1])

print("ANOVA method:",
      anova_matrix.shape[1])


# =================================================
# 2. COMMON FEATURES
# =================================================

variance_features = set(variance_matrix.columns)
anova_features = set(anova_matrix.columns)

common_features = variance_features.intersection(
    anova_features
)

print("\nCommon features between methods:",
      len(common_features))


# =================================================
# 3. UNIQUE FEATURES
# =================================================

variance_unique = variance_features - anova_features
anova_unique = anova_features - variance_features

print("Variance-only features:",
      len(variance_unique))

print("ANOVA-only features:",
      len(anova_unique))


# =================================================
# 4. PCA EXPLAINED VARIANCE COMPARISON
# =================================================

variance_pca_score = (
    pca_variance.explained_variance_ratio_.sum()
)

anova_pca_score = (
    pca_anova.explained_variance_ratio_.sum()
)

print("\nPCA explained variance:")

print(
    "Variance method:",
    variance_pca_score
)

print(
    "ANOVA method:",
    anova_pca_score
)


# =================================================
# 5. COMPARISON
# =================================================

print("\n================================")
print("FINAL COMPARISON")
print("================================")

if anova_pca_score > variance_pca_score:

    print(
        "ANOVA feature selection provides "
        "higher PCA explained variance."
    )

    print(
        "Therefore, ANOVA performed better "
        "for representing group-related variation."
    )

elif variance_pca_score > anova_pca_score:

    print(
        "Variance feature selection provides "
        "higher PCA explained variance."
    )

    print(
        "Therefore, Variance selection performed "
        "better for representing variation."
    )

else:

    print(
        "Both methods provide the same "
        "PCA explained variance."
    )


# =================================================
# 6. FEATURE OVERLAP
# =================================================

total_selected = 5000

overlap_percentage = (
    len(common_features) /
    total_selected
) * 100

print("\nFeature overlap:")
print(
    "Common features:",
    len(common_features)
)

print(
    "Overlap percentage:",
    overlap_percentage,
    "%"
)


print("\n================================")
print("STEP 9 COMPLETED")
print("================================")
