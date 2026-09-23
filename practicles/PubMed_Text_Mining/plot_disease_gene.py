import pandas as pd
import matplotlib.pyplot as plt

# Load disease-gene results
df = pd.read_csv("disease_gene_signals.csv")

# Create label for plotting
df["Association"] = df["Disease"] + " - " + df["Gene"]

# Select top 10 associations
top10 = df.head(10)

# Plot
plt.figure(figsize=(10, 6))

plt.barh(
    top10["Association"],
    top10["Co_occurrence"]
)

plt.xlabel("Co-occurrence Count")
plt.ylabel("Disease - Gene")
plt.title("Top Disease-Gene Associations from PubMed Abstracts")

plt.gca().invert_yaxis()

plt.tight_layout()

# Save figure
plt.savefig("disease_gene_signals.png", dpi=300)

plt.show()

print("Plot saved as: disease_gene_signals.png")

