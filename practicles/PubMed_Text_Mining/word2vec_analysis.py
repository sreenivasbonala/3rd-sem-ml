import pandas as pd
import re

from gensim.models import Word2Vec

# 1. Load PubMed abstracts
df = pd.read_csv("pubmed_abstracts.csv")

# 2. Remove missing abstracts
df = df.dropna(subset=["Abstract"])

# 3. Convert abstracts into sentences of words
sentences = []

for abstract in df["Abstract"]:
    abstract = abstract.lower()
    abstract = re.sub(r"[^a-zA-Z\s]", " ", abstract)

    words = abstract.split()

    # Remove very short words
    words = [word for word in words if len(word) > 2]

    sentences.append(words)

print("Number of abstracts:", len(sentences))
print("First abstract words:")
print(sentences[0][:20])

# 4. Train Word2Vec model
model = Word2Vec(
    sentences=sentences,
    vector_size=100,
    window=5,
    min_count=2,
    workers=4,
    sg=1
)

# 5. Save the model
model.save("pubmed_word2vec.model")

print("\nWord2Vec model trained successfully!")
print("Vocabulary size:", len(model.wv))
print("Model saved as: pubmed_word2vec.model")

