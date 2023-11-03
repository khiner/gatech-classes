import cvxpy as cp
import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split

data = load_breast_cancer()
X = 2 * data.data - 1
y = 2 * data.target - 1
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

N, D = X_train.shape
c = 1 / 3
r = 1

w = cp.Variable(D)
b = cp.Variable()
t = cp.Variable(N)

# Set up the constraints
constraints = []
for i in range(N):
    score = y_train[i] * (X_train[i, :] @ w + b)
    constraints.append(t[i] >= 1 - ((1 - c) * score / c))
    constraints.append(t[i] >= 1 - score)
    constraints.append(t[i] >= 0)

constraints.append(cp.norm(w) <= r)

# Set up the objective
objective = cp.Minimize(cp.sum(t) / N)

# Formulate and solve the problem
problem = cp.Problem(objective, constraints)
problem.solve()

# Extract the optimal values
w_opt = w.value
b_opt = b.value

# print(f"Optimal w: {w_opt}")
# print(f"Optimal b: {b_opt}")

def custom_loss(y, scores, c):
    loss_values = np.zeros(scores.shape)
    # Condition: y * score < 0
    mask1 = y * scores < 0
    loss_values[mask1] = 1 - ((1 - c) * y[mask1] * scores[mask1]) / c
    # Condition: 0 <= y * score <= 1
    mask2 = (y * scores >= 0) & (y * scores <= 1)
    loss_values[mask2] = 1 - y[mask2] * scores[mask2]
    # Condition: y * score > 1
    loss_values[y * scores > 1] = 0
    
    return np.mean(loss_values)

scores_test = X_test @ w_opt + b_opt
test_loss = custom_loss(y_test, scores_test, c)

print(f"Test loss: {test_loss}")
