import pandas as pd
import itertools
import matplotlib.pyplot as plt

# Load expression data
df = pd.read_csv("expression_matrix.csv")

# Remove probe ID column
expression = df.drop(columns=["ID_REF"])

# Select 15 most variable probes
variances = expression.var(axis=1)
top_probes = variances.nlargest(15).index

selected = expression.loc[top_probes]

# Transpose: samples x probes
selected = selected.T

# Convert expression values to High/Low using median
binary = selected.apply(
    lambda col: col >= col.median()
).astype(int)

print("Binary matrix:", binary.shape)

# ------------------------------------------------
# Pairwise association rules
# ------------------------------------------------

rules = []

for A, B in itertools.permutations(binary.columns, 2):

    # A and B are both high
    both = ((binary[A] == 1) & (binary[B] == 1)).sum()

    # A is high
    A_count = (binary[A] == 1).sum()

    # B is high
    B_count = (binary[B] == 1).sum()

    total = len(binary)

    support = both / total
    confidence = both / A_count

    # Lift
    expected = (A_count / total) * (B_count / total)

    if expected > 0:
        lift = support / expected
    else:
        lift = 0

    # Keep support >= 40% and confidence >= 70%
    if support >= 0.40 and confidence >= 0.70:
        rules.append([
            A,
            B,
            support,
            confidence,
            lift
        ])

# Create dataframe
rules_df = pd.DataFrame(
    rules,
    columns=[
        "Antecedent",
        "Consequent",
        "Support",
        "Confidence",
        "Lift"
    ]
)

# Sort by lift
rules_df = rules_df.sort_values(
    by="Lift",
    ascending=False
)

# Save results
rules_df.to_csv(
    "association_rules.csv",
    index=False
)

print("\nAssociation Rules:")
print(rules_df.head(20))

print("\nTotal rules:", len(rules_df))

# ------------------------------------------------
# Plot strongest rules
# ------------------------------------------------

if len(rules_df) > 0:

    top = rules_df.head(10).copy()

    top["Rule"] = (
        top["Antecedent"].astype(str)

    + " -> "

    + top["Consequent"].astype(str)

)

    plt.figure(figsize=(10, 6))

    plt.barh(
        top["Rule"],
        top["Lift"]
    )

    plt.xlabel("Lift")
    plt.ylabel("Gene/Probe Association")
    plt.title("Top Association Rules")

    plt.gca().invert_yaxis()

    plt.tight_layout()

    plt.savefig(
        "association_rules.png",
        dpi=300
    )

    print("\nPlot saved as association_rules.png")

print("\nAnalysis completed successfully.")

