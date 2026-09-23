# -*- coding: utf-8 -*-
"""
Basics of ML - Iris Classification
Log Transformation + Before/After Plots + Logistic Regression
"""

# ==========================================
# STEP 1: IMPORT REQUIRED LIBRARIES
# ==========================================

from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ==========================================
# STEP 2: LOAD AND EXPLORE THE IRIS DATA
# ==========================================

iris = load_iris()

# Convert dataset into Pandas DataFrame
df = pd.DataFrame(
    data=iris.data,
    columns=iris.feature_names
)

df["target_flower"] = iris.target

print("\nFirst 5 rows of the dataset:")
print(df.head())

print("\nFlower types to predict:")
print(list(iris.target_names))


# ==========================================
# STEP 3: EXTRACT FEATURES AND TARGET
# ==========================================

X = iris.data
y = iris.target

print("\nOriginal X shape:")
print(X.shape)

print("\nOriginal y shape:")
print(y.shape)


# ==========================================
# STEP 4: LOG TRANSFORMATION
# ==========================================

# Keep original data
X_original = X.copy()

# Apply Log2 transformation
X_log = np.log2(X_original + 1)

print("\nOriginal Data:")
print(X_original[:5])

print("\nAfter Log Transformation:")
print(X_log[:5])


# ==========================================
# STEP 5: PLOT SEPAL LENGTH
# BEFORE LOG TRANSFORMATION
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(X_original[:, 0], bins=15)

plt.xlabel("Sepal Length")
plt.ylabel("Frequency")
plt.title("Sepal Length Before Log Transformation")

plt.show()


# ==========================================
# STEP 6: PLOT SEPAL LENGTH
# AFTER LOG TRANSFORMATION
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(X_log[:, 0], bins=15)

plt.xlabel("log2(Sepal Length + 1)")
plt.ylabel("Frequency")
plt.title("Sepal Length After Log Transformation")

plt.show()


# ==========================================
# STEP 7: PLOT PETAL LENGTH
# BEFORE LOG TRANSFORMATION
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(X_original[:, 2], bins=15)

plt.xlabel("Petal Length")
plt.ylabel("Frequency")
plt.title("Petal Length Before Log Transformation")

plt.show()


# ==========================================
# STEP 8: PLOT PETAL LENGTH
# AFTER LOG TRANSFORMATION
# ==========================================

plt.figure(figsize=(8, 5))

plt.hist(X_log[:, 2], bins=15)

plt.xlabel("log2(Petal Length + 1)")
plt.ylabel("Frequency")
plt.title("Petal Length After Log Transformation")

plt.show()


# ==========================================
# STEP 9: USE LOG-TRANSFORMED DATA FOR ML
# ==========================================

X = X_log

print("\nLog-transformed X shape:")
print(X.shape)


# ==========================================
# STEP 10: TRAIN-TEST SPLIT
# ==========================================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.50,
    random_state=1
)

print("\nData Split Complete!")

print("Training samples:", X_train.shape[0])
print("Testing samples:", X_test.shape[0])


# ==========================================
# STEP 11: INITIALIZE LOGISTIC REGRESSION
# ==========================================

model = LogisticRegression(
    max_iter=200
)

print("\nTraining Logistic Regression model...")

model.fit(X_train, y_train)

print("Training complete!")


# ==========================================
# STEP 12: MAKE PREDICTIONS
# ==========================================

print("\nTesting the model on unseen data...")

predictions = model.predict(X_test)


# ==========================================
# STEP 13: CALCULATE ACCURACY
# ==========================================

accuracy = accuracy_score(
    y_test,
    predictions
)

print("\nFinal Accuracy Score:")
print(f"{accuracy * 100:.2f}%")


# ==========================================
# STEP 14: CLASSIFICATION REPORT
# ==========================================

print("\nDetailed Performance Report:")

print(
    classification_report(
        y_test,
        predictions,
        target_names=iris.target_names
    )
)


# ==========================================
# STEP 15: CUSTOM FLOWER PREDICTION
# ==========================================

print("\nPredicting a custom flower measurement...")

# Original measurements
my_custom_flower = [
    [5.1, 3.5, 1.4, 0.2]
]

# Apply the SAME log transformation
my_custom_flower_log = np.log2(
    np.array(my_custom_flower) + 1
)

# Prediction
predicted_class = model.predict(
    my_custom_flower_log
)

# Convert class number to flower name
predicted_name = iris.target_names[
    predicted_class[0]
]


# ==========================================
# STEP 16: DISPLAY CUSTOM PREDICTION
# ==========================================

print("\nInput measurements:")
print(my_custom_flower[0])

print("\nLog-transformed input:")
print(my_custom_flower_log[0])

print(
    f"\nThe model predicts this flower is: "
    f"{predicted_name.upper()}"
)


# ==========================================
# STEP 17: PREDICTION PROBABILITY
# ==========================================

probabilities = model.predict_proba(
    my_custom_flower_log
)

print("\nPrediction probabilities:")

for flower, probability in zip(
    iris.target_names,
    probabilities[0]
):
    print(
        f"{flower}: "
        f"{probability * 100:.2f}%"
    )


print("\nProgram completed successfully!")