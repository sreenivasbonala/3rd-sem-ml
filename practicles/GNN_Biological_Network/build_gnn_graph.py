import pandas as pd
import torch
from torch_geometric.data import Data
from sklearn.preprocessing import StandardScaler

# -----------------------------
# Load gene-disease information
# -----------------------------
gene_df = pd.read_csv("gene_disease_clean.csv")

# Load STRING interactions
edge_df = pd.read_csv("string_edges.csv")

# All 353 genes
genes = gene_df["Gene"].dropna().drop_duplicates().tolist()

gene_to_id = {gene: i for i, gene in enumerate(genes)}

# -----------------------------
# Create node features
# -----------------------------
feature_cols = [
    "score",
    "Normalized score",
    "N PMIDs",
    "N CTs",
    "N variants gda",
    "DSI g",
    "DPI g",
    "pLI"
]

features = gene_df[feature_cols].apply(
    pd.to_numeric, errors="coerce"
).fillna(0)

# Standardize features
scaler = StandardScaler()
X = scaler.fit_transform(features)

x = torch.tensor(X, dtype=torch.float)

# -----------------------------
# Create graph edges
# -----------------------------
edges = []

for _, row in edge_df.iterrows():

    gene_a = row["preferredName_A"]
    gene_b = row["preferredName_B"]

    if gene_a in gene_to_id and gene_b in gene_to_id:

        a = gene_to_id[gene_a]
        b = gene_to_id[gene_b]

        edges.append([a, b])
        edges.append([b, a])

edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()

# -----------------------------
# Create PyTorch Geometric Data
# -----------------------------
data = Data(
    x=x,
    edge_index=edge_index
)

print("\n==============================")
print("GNN GRAPH CREATED")
print("==============================")
print("Number of nodes:", data.num_nodes)
print("Number of edges:", data.num_edges)
print("Node feature dimensions:", data.num_node_features)
print("Feature matrix shape:", data.x.shape)
print("Edge index shape:", data.edge_index.shape)

print("\nGenes participating in STRING network:",
      len(set(edge_df["preferredName_A"]) |
          set(edge_df["preferredName_B"])))

print("Genes without STRING edges:",
      len(genes) -
      len(set(edge_df["preferredName_A"]) |
          set(edge_df["preferredName_B"])))

# Save tensors
torch.save(data, "gene_network.pt")

# Save gene mapping
pd.DataFrame({
    "node_id": range(len(genes)),
    "Gene": genes
}).to_csv("gene_node_mapping.csv", index=False)

print("\nSaved:")
print("gene_network.pt")
print("gene_node_mapping.csv")
