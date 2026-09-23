from Bio import Entrez, Medline
import pandas as pd

# Tell NCBI who is making the request
Entrez.email = "bonala980@gmail.com"

# PubMed search query
query = 'cancer AND gene association'

# Search PubMed
handle = Entrez.esearch(
    db="pubmed",
    term=query,
    retmax=100
)

record = Entrez.read(handle)
handle.close()

pmids = record["IdList"]

print("Number of papers found:", len(pmids))

# Download abstracts
handle = Entrez.efetch(
    db="pubmed",
    id=pmids,
    rettype="medline",
    retmode="text"
)

records = list(Medline.parse(handle))
handle.close()

# Extract useful information
data = []

for record in records:
    data.append({
        "PMID": record.get("PMID", ""),
        "Title": record.get("TI", ""),
        "Abstract": record.get("AB", ""),
        "Authors": record.get("AU", ""),
        "Journal": record.get("JT", "")
    })

# Save as CSV
df = pd.DataFrame(data)

df.to_csv("pubmed_abstracts.csv", index=False)

print("Saved:", len(df), "records")
print("File: pubmed_abstracts.csv")

