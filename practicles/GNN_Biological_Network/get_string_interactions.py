import pandas as pd
import requests

# Load the 353 genes
genes = pd.read_csv("genes.txt", header=None)[0].dropna().tolist()

print("Genes submitted to STRING:", len(genes))

# STRING API
url = "https://string-db.org/api/tsv/network"

params = {
    "identifiers": "%0d".join(genes),
    "species": 9606,
    "required_score": 700,
    "caller_identity": "GNN_Biological_Network"
}

response = requests.get(url, params=params, timeout=120)

print("HTTP status:", response.status_code)

if response.status_code != 200:
    print("STRING request failed")
    print(response.text[:500])
    raise SystemExit(1)

with open("string_interactions.tsv", "w") as f:
    f.write(response.text)

print("STRING interaction file saved: string_interactions.tsv")

# Read and inspect
df = pd.read_csv("string_interactions.tsv", sep="\t")

print("Interactions retrieved:", len(df))
print("Columns:")
print(df.columns.tolist())

print("\nFirst 5 interactions:")
print(df.head().to_string(index=False))
