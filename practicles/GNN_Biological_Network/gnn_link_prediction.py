
import torch
import pandas as pd
import matplotlib.pyplot as plt

from torch import nn
from torch_geometric.nn import GCNConv
from torch_geometric.transforms import RandomLinkSplit
from sklearn.metrics import roc_auc_score, average_precision_score


# ============================================================
# 1. Load biological graph
# ============================================================

data = torch.load("gene_network.pt", weights_only=False)

print("\n====================================")
print("GCN LINK PREDICTION")
print("====================================")
print("Nodes:", data.num_nodes)
print("Edges:", data.num_edges)
print("Node features:", data.num_node_features)


# ============================================================
# 2. Train / validation / test split
# ============================================================

transform = RandomLinkSplit(
    num_val=0.10,
    num_test=0.20,
    is_undirected=True,
    add_negative_train_samples=True
)

train_data, val_data, test_data = transform(data)

print("\nGraph split:")
print("Training edges:", train_data.edge_index.size(1))
print("Validation edges:", val_data.edge_label_index.size(1))
print("Test edges:", test_data.edge_label_index.size(1))


# ============================================================
# 3. GCN model
# ============================================================

class GCNLinkPredictor(nn.Module):

    def __init__(self, input_dim, hidden_dim=32, embedding_dim=16):

        super().__init__()

        self.conv1 = GCNConv(input_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, embedding_dim)

    def encode(self, x, edge_index):

        x = self.conv1(x, edge_index)
        x = torch.relu(x)

        x = self.conv2(x, edge_index)

        return x

    def decode(self, z, edge_index):

        source = z[edge_index[0]]
        target = z[edge_index[1]]

        return (source * target).sum(dim=1)


# ============================================================
# 4. Device
# ============================================================

if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("\nUsing Apple Silicon MPS")
else:
    device = torch.device("cpu")
    print("\nUsing CPU")


train_data = train_data.to(device)
val_data = val_data.to(device)
test_data = test_data.to(device)


# ============================================================
# 5. Initialize model
# ============================================================

model = GCNLinkPredictor(
    input_dim=data.num_node_features
).to(device)

optimizer = torch.optim.Adam(
    model.parameters(),
    lr=0.01,
    weight_decay=5e-4
)

criterion = nn.BCEWithLogitsLoss()


# ============================================================
# 6. Training
# ============================================================

losses = []

epochs = 200

for epoch in range(1, epochs + 1):

    model.train()

    optimizer.zero_grad()

    z = model.encode(
        train_data.x,
        train_data.edge_index
    )

    pred = model.decode(
        z,
        train_data.edge_label_index
    )

    loss = criterion(
        pred,
        train_data.edge_label.float()
    )

    loss.backward()

    optimizer.step()

    losses.append(loss.item())

    if epoch == 1 or epoch % 20 == 0:

        print(
            f"Epoch {epoch:03d} | "
            f"Loss: {loss.item():.4f}"
        )


# ============================================================
# 7. Evaluation function
# ============================================================

def evaluate(model, graph):

    model.eval()

    with torch.no_grad():

        z = model.encode(
            graph.x,
            graph.edge_index
        )

        logits = model.decode(
            z,
            graph.edge_label_index
        )

        probabilities = torch.sigmoid(logits)

        y_true = graph.edge_label.cpu().numpy()

        y_score = probabilities.cpu().numpy()

        auc = roc_auc_score(
            y_true,
            y_score
        )

        ap = average_precision_score(
            y_true,
            y_score
        )

    return auc, ap


# ============================================================
# 8. Validation and test
# ============================================================

val_auc, val_ap = evaluate(
    model,
    val_data
)

test_auc, test_ap = evaluate(
    model,
    test_data
)


# ============================================================
# 9. Results
# ============================================================

print("\n====================================")
print("RESULTS")
print("====================================")

print(f"Validation ROC-AUC: {val_auc:.4f}")
print(f"Validation Average Precision: {val_ap:.4f}")

print(f"Test ROC-AUC:       {test_auc:.4f}")
print(f"Test Average Precision: {test_ap:.4f}")


# ============================================================
# 10. Save training curve
# ============================================================

plt.figure(figsize=(8, 5))

plt.plot(
    range(1, epochs + 1),
    losses
)

plt.xlabel("Epoch")
plt.ylabel("Training Loss")
plt.title("GCN Link Prediction Training Loss")

plt.tight_layout()

plt.savefig(
    "gnn_training_curve.png",
    dpi=300
)

plt.close()


# ============================================================
# 11. Save results
# ============================================================

results = pd.DataFrame({
    "Metric": [
        "Nodes",
        "Original directed edges",
        "Training edges",
        "Validation edges",
        "Test edges",
        "Node features",
        "Validation ROC-AUC",
        "Validation Average Precision",
        "Test ROC-AUC",
        "Test Average Precision"
    ],

    "Value": [
        data.num_nodes,
        data.num_edges,
        train_data.edge_index.size(1),
        val_data.edge_label_index.size(1),
        test_data.edge_label_index.size(1),
        data.num_node_features,
        val_auc,
        val_ap,
        test_auc,
        test_ap
    ]
})

results.to_csv(
    "gnn_link_prediction_results.csv",
    index=False
)

print("\nSaved:")
print("gnn_training_curve.png")
print("gnn_link_prediction_results.csv")

print("\n====================================")
print("ANALYSIS COMPLETED")
print("====================================")

