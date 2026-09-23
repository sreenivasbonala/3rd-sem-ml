import pandas as pd
import re
from collections import Counter

# Load PubMed data
df = pd.read_csv("pubmed_abstracts.csv")

# Remove abstracts that are missing
df = df.dropna(subset=["Abstract"])

print("Number of abstracts:", len(df))

# Disease terms
diseases = [
    "cancer",
    "breast cancer",
    "lung cancer",
    "prostate cancer",
    "colorectal cancer",
    "melanoma",
    "liver cancer",
    "pancreatic cancer",
    "ovarian cancer"
]

# Important cancer-related gene symbols
genes = [
    "TP53",
    "BRCA1",
    "BRCA2",
    "EGFR",
    "KRAS",
    "MYC",
    "VEGFA",
    "METTL3",
    "IGF2BP2",
    "PTEN",
    "PIK3CA",
    "AKT1",
    "ERBB2",
    "BRAF",
    "ALK"
]

# Store disease-gene associations
results = []

# Analyze every abstract
for abstract in df["Abstract"]:

    abstract_text = str(abstract)

    # Search for diseases
    found_diseases = []

    for disease in diseases:
        if re.search(r"\b" + re.escape(disease) + r"\b",
                     abstract_text,
                     re.IGNORECASE):
            found_diseases.append(disease)

    # Search for genes
    found_genes = []

    for gene in genes:
        if re.search(r"\b" + re.escape(gene) + r"\b",
                     abstract_text):
            found_genes.append(gene)

    # Create disease-gene pairs
    for disease in found_diseases:
        for gene in found_genes:
            results.append((disease, gene))


# Count how many abstracts contain each pair
pair_counts = Counter(results)

# Convert results into a DataFrame
output = []

for (disease, gene), count in pair_counts.items():
    output.append({
        "Disease": disease,
        "Gene": gene,
        "Co_occurrence": count
    })

result_df = pd.DataFrame(output)

# Sort by strongest association
if not result_df.empty:
    result_df = result_df.sort_values(
        by="Co_occurrence",
        ascending=False
    )

# Save results
result_df.to_csv(
    "disease_gene_signals.csv",
    index=False
)

# Display top results
print("\nTop Disease-Gene Signals:")
print(result_df.head(20).to_string(index=False))

print("\nSaved: disease_gene_signals.csv")
