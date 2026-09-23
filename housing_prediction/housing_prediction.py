import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

# 1. Synthetic historical data
np.random.seed(42)

years = np.arange(2016, 2026)

prices = np.array([
    45, 48, 52, 56, 60,
    65, 71, 77, 83, 90
])

# 2. Create DataFrame
data = pd.DataFrame({
    "Year": years,
    "Price_Lakhs": prices
})

print(data)

# 3. Feature and target
X = data[["Year"]]
y = data["Price_Lakhs"]

# 4. Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42
)

# 5. Create model
model = LinearRegression()

# 6. Train model
model.fit(X_train, y_train)

# 7. Test prediction
y_pred = model.predict(X_test)

# 8. Calculate MSE
mse = mean_squared_error(y_test, y_pred)

print("\nSlope:", model.coef_[0])
print("Intercept:", model.intercept_)
print("MSE:", mse)

# 9. Future years
future_years = np.arange(2026, 2036).reshape(-1, 1)

future_prices = model.predict(future_years)

# 10. Display future predictions
future_data = pd.DataFrame({
    "Year": future_years.flatten(),
    "Predicted_Price_Lakhs": future_prices
})

print("\nFuture Housing Price Prediction:")
print(future_data)

# 11. Plot
plt.scatter(years, prices, label="Historical Data")
plt.plot(years, model.predict(X), label="Regression Line")
plt.plot(
    future_years,
    future_prices,
    marker="o",
    label="Future Prediction"
)

plt.xlabel("Year")
plt.ylabel("House Price (₹ Lakhs)")
plt.title("Housing Price Prediction")
plt.legend()
plt.show()