import pandas as pd
import re
import nltk

from nltk.corpus import stopwords
from sklearn.feature_extraction.text import TfidfVectorizer

# Download stopwords
nltk.download("stopwords")

# Load PubMed data
df = pd.read_csv("pubmed_abstracts.csv")

# Remove missing abstracts
df = df.dropna(subset=["Abstract"])

# Text preprocessing
stop_words = set(stopwords.words("english"))

def preprocess(text):
    text = text.lower()
    text = re.sub(r"[^a-zA-Z\s]", " ", text)
    words = text.split()
    words = [word for word in words if word not in stop_words]
    return " ".join(words)

df["Clean_Abstract"] = df["Abstract"].apply(preprocess)

# TF-IDF
vectorizer = TfidfVectorizer(
    max_features=100
)

tfidf_matrix = vectorizer.fit_transform(df["Clean_Abstract"])

# Get terms
terms = vectorizer.get_feature_names_out()

# Calculate average TF-IDF score
scores = tfidf_matrix.mean(axis=0).A1

results = pd.DataFrame({
    "Term": terms,
    "TF_IDF_Score": scores
})

# Sort by importance
results = results.sort_values(
    "TF_IDF_Score",
    ascending=False
)

# Save results
results.to_csv("tfidf_results.csv", index=False)

print("\nTop 20 important terms:")
print(results.head(20))

print("\nSaved: tfidf_results.csv")
