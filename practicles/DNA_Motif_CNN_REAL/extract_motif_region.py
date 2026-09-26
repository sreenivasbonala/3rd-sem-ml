import pandas as pd

# Load test predictions
df = pd.read_csv("real_test_predictions.csv")

# Select correctly predicted positive CTCF sequences
positive = df[
    (df["label"] == 1) &
    (df["predicted_label"] == 1)
].copy()

print("Correctly predicted positive sequences:", len(positive))

# Extract region around the important CNN positions
start = 90
end = 121

print("\n===== CTCF IMPORTANT REGION =====")
print(f"Positions: {start+1}-{end}")

for i, row in positive.head(10).iterrows():

    seq = row["sequence"]

    region = seq[start:end]

    print(f"\n{row['id']}")
    print(f"Full length : {len(seq)}")
    print(f"Region      : {region}")

print("\nDone.")
