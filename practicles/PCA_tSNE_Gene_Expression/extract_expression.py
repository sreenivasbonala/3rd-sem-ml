import GEOparse
import pandas as pd

print("Loading GSE10072...")

g = GEOparse.get_GEO("GSE10072", destdir=".")

samples = []

for gsm_id, gsm in g.gsms.items():
    table = gsm.table
    table = table.set_index("ID_REF")
    samples.append(table["VALUE"].rename(gsm_id))

expression = pd.concat(samples, axis=1)

print("Expression matrix shape:", expression.shape)

expression.to_csv("expression_matrix.csv")

print("Saved: expression_matrix.csv")
