import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# ============================================================
# 1. LOAD GENE-DISEASE DATA
# ============================================================

df = pd.read_csv("gene_disease_clean.csv")

print("Dataset shape:", df.shape)

print("\nFirst 5 rows:")
print(df.head())


# ============================================================
# 2. CREATE HETEROGENEOUS GENE-DISEASE NETWORK
# ============================================================

G = nx.Graph()

disease = "Stomach Neoplasms"

# Add disease node
G.add_node(
    disease,
    node_type="Disease"
)

# Add genes and gene-disease associations
for _, row in df.iterrows():

    gene = str(row["Gene"])

    # Add gene node
    G.add_node(
        gene,
        node_type="Gene"
    )

    # Add gene-disease edge
    G.add_edge(
        gene,
        disease,
        edge_type="GDA",
        score=float(row["score"])
    )


# ============================================================
# 3. NETWORK INFORMATION
# ============================================================

print("\nNetwork Information")
print("-------------------")

print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())

gene_nodes = [
    node for node, data in G.nodes(data=True)
    if data["node_type"] == "Gene"
]

disease_nodes = [
    node for node, data in G.nodes(data=True)
    if data["node_type"] == "Disease"
]

print("Gene nodes:", len(gene_nodes))
print("Disease nodes:", len(disease_nodes))


# ============================================================
# 4. CENTRALITY ANALYSIS
# ============================================================

degree = nx.degree_centrality(G)

betweenness = nx.betweenness_centrality(G)

closeness = nx.closeness_centrality(G)


# Create results table
results = []

for gene in gene_nodes:

    results.append({
        "Gene": gene,
        "Degree_Centrality": degree[gene],
        "Betweenness_Centrality": betweenness[gene],
        "Closeness_Centrality": closeness[gene]
    })


centrality_df = pd.DataFrame(results)


# Sort by degree centrality
centrality_df = centrality_df.sort_values(
    "Degree_Centrality",
    ascending=False
)


# Save results
centrality_df.to_csv(
    "centrality_results.csv",
    index=False
)


print("\nTop 10 genes by centrality:")
print(
    centrality_df.head(10).to_string(index=False)
)


# ============================================================
# 5. COMMUNITY DETECTION
# ============================================================

communities = nx.community.greedy_modularity_communities(G)

print("\nNumber of communities:", len(communities))

for i, community in enumerate(communities, start=1):

    print(
        f"Community {i}:",
        list(community)[:10]
    )


# ============================================================
# 6. NETWORK VISUALIZATION
# ============================================================

print("\nCreating network visualization...")

plt.figure(figsize=(16, 14))


# IMPORTANT:
# Generate positions for ALL nodes
pos = nx.spring_layout(
    G,
    seed=42
)


# Draw gene nodes
nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=gene_nodes,
    node_size=250,
    node_color="skyblue",
    node_shape="o"
)


# Draw disease node
nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=disease_nodes,
    node_size=1800,
    node_color="red",
    node_shape="s"
)


# Draw edges
nx.draw_networkx_edges(
    G,
    pos,
    alpha=0.4
)


# Draw labels
nx.draw_networkx_labels(
    G,
    pos,
    font_size=6
)


plt.title(
    "Gene-Disease Network: Stomach Neoplasms",
    fontsize=16
)

plt.axis("off")

plt.tight_layout()


# Save figure
plt.savefig(
    "gene_disease_network.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 7. FINAL MESSAGE
# ============================================================

print("\n====================================")
print("ANALYSIS COMPLETED SUCCESSFULLY")
print("====================================")

print("Network plot:")
print("gene_disease_network.png")

print("\nCentrality results:")
print("centrality_results.csv")
