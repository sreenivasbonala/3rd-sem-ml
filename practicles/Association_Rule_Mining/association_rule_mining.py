import pandas as pd
import matplotlib.pyplot as plt
from mlxtend.frequent_patterns import apriori, association_rules

# Load expression matrix
df = pd.read_csv("expression_matrix.csv")

# Separate gene IDs and expression values
gene_ids = df["ID_REF"]
expression = df.drop(columns=["ID_REF"])

print("Original expression matrix:", expression.shape)

# Select the 30 most variable genes
variances = expression.var(axis=1)
top_indices = variances.nlargest(30).index

selected = expression.loc[top_indices].copy()
selected["ID_REF"] = gene_ids.loc[top_indices].values

# Transpose: samples × genes
sample_gene = selected.set_index("ID_REF").T

print("Selected matrix:", sample_gene.shape)

# Convert continuous expression into binary:
# 1 = High expression
# 0 = Low expression
binary = sample_gene.apply(
    lambda col: col >= col.median()
).astype(bool)

print("Binary matrix:", binary.shape)

# Apriori frequent itemsets
frequent_itemsets = apriori(
    binary,
    min_support=0.30,
    use_colnames=True
)

frequent_itemsets["itemsets"] = frequent_itemsets["itemsets"].apply(
    lambda x: ", ".join(sorted(x))
)

frequent_itemsets.to_csv(
    "frequent_itemsets.csv",
    index=False
)

print("\nFrequent itemsets:")
print(frequent_itemsets.head(10))

# Association rules
# Re-create itemsets with frozensets for mlxtend
freq_for_rules = apriori(
    binary,
    min_support=0.30,
    use_colnames=True
)

rules = association_rules(
    freq_for_rules,
    metric="confidence",
    min_threshold=0.70
)

if not rules.empty:
    rules["antecedents"] = rules["antecedents"].apply(
        lambda x: ", ".join(sorted(x))
    )
    rules["consequents"] = rules["consequents"].apply(
        lambda x: ", ".join(sorted(x))
    )

    rules = rules.sort_values(
        ["lift", "confidence"],
        ascending=False
    )

    rules.to_csv(
        "association_rules.csv",
        index=False
    )

    print("\nTop association rules:")
    print(
        rules[
            [
                "antecedents",
                "consequents",
                "support",
                "confidence",
                "lift"
            ]
        ].head(10)
    )

    # Plot top 10 rules by lift
    top_rules = rules.head(10).copy()
    labels = (
        top_rules["antecedents"]
        + " → "
        + top_rules["consequents"]
    )

    plt.figure(figsize=(10, 6))
    plt.barh(labels, top_rules["lift"])
    plt.xlabel("Lift")
    plt.ylabel("Association Rule")
    plt.title("Top Gene Association Rules")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("association_rules.png", dpi=300)
    plt.close()

    print("\nSaved:")
    print("association_rules.csv")
    print("association_rules.png")

else:
    print("\nNo association rules found.")
    print("Try lowering min_support or min_threshold.")
