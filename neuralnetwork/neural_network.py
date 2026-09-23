import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

from sklearn.datasets import make_moons
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


# ==========================================
# 1. CREATE DATASET
# ==========================================

X, y = make_moons(
    n_samples=1000,
    noise=0.2,
    random_state=42
)

print("Original X shape:", X.shape)
print("Original y shape:", y.shape)


# ==========================================
# 2. TRAIN-TEST SPLIT
# ==========================================

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42
)

print("Training samples:", X_train.shape[0])
print("Testing samples:", X_test.shape[0])


# ==========================================
# 3. STANDARDIZATION
# ==========================================

scaler = StandardScaler()

X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)


# ==========================================
# 4. BUILD NEURAL NETWORK
# ==========================================

model = Sequential([
    
    Dense(
        16,
        input_dim=2,
        activation="relu"
    ),

    Dense(
        8,
        activation="relu"
    ),

    Dense(
        1,
        activation="sigmoid"
    )
])


# ==========================================
# 5. COMPILE MODEL
# ==========================================

model.compile(
    optimizer="adam",
    loss="binary_crossentropy",
    metrics=["accuracy"]
)


# ==========================================
# 6. SHOW MODEL STRUCTURE
# ==========================================

model.summary()


# ==========================================
# 7. TRAIN MODEL
# ==========================================

history = model.fit(
    X_train,
    y_train,
    epochs=50,
    batch_size=32,
    validation_data=(X_test, y_test)
)


# ==========================================
# 8. EVALUATE MODEL
# ==========================================

loss, accuracy = model.evaluate(
    X_test,
    y_test
)

print(
    f"\nTest Accuracy: {accuracy * 100:.2f}%"
)