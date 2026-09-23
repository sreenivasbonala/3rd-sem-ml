import GEOparse
import pandas as pd

print("Loading GSE10072...")

g = GEOparse.get_GEO("GSE10072", destdir=".")

rows = []

for gsm_id, gsm in g.gsms.items():

    title = gsm.metadata.get("title", [""])[0]

    if "Lung Tumor" in title:
        group = "Lung Tumor"
    elif "Normal Lung" in title:
        group = "Normal Lung"
    else:
        group = "Unknown"

    rows.append({
        "Sample": gsm_id,
        "Title": title,
        "Group": group
    })

metadata = pd.DataFrame(rows)

metadata.to_csv("sample_metadata.csv", index=False)

print("\nGroup counts:")
print(metadata["Group"].value_counts())

print("\nSaved: sample_metadata.csv")

