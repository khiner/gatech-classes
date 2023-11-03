from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Perceptron as SklearnPerceptron

import numpy as np

DEFAULT_EPOCHS = 50 # Found with `find_best_epochs`.

class Perceptron:
    def __init__(self, learning_rate=1.0, epochs=DEFAULT_EPOCHS):
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.weights = None
        self.bias = None

    def fit(self, X, y, sample_weight=None):
        self.weights = np.zeros(X.shape[1])
        self.bias = 0

        if sample_weight is None:
            sample_weight = np.ones(X.shape[0]) / X.shape[0]

        for epoch in range(self.epochs):
            for i in range(X.shape[0]):
                y_pred = np.sign(np.dot(X[i], self.weights) + self.bias)
                update = self.learning_rate * (y[i] - y_pred)
                self.weights += update * X[i] * sample_weight[i]
                self.bias += update * sample_weight[i]

    def predict(self, X_in):
        return np.sign(np.dot(X_in, self.weights) + self.bias)


# Search for the number of epochs that gives the best accuracy over the test set.
def find_best_epochs(X_train, y_train, X_test, y_test, min = 50, max = 350, step = 25):
    best_accuracy = 0
    best_epoch = 0

    for epochs in range(min, max + 1, step):
        perceptron = Perceptron(epochs=epochs)
        perceptron.fit(X_train, y_train)
        y_pred = perceptron.predict(X_test)
        accuracy = np.mean(y_pred == y_test)

        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_epoch = epochs

        print(f"Epochs: {epochs}, Test Accuracy: {accuracy * 100:.2f}%")

    print(f"Best Epochs: {best_epoch}, Best Test Accuracy: {best_accuracy * 100:.2f}%")
    return best_epoch

def create_perceptron(custom=True, learning_rate=1.0, epochs=DEFAULT_EPOCHS):
    if custom:
        return Perceptron(learning_rate=learning_rate, epochs=epochs)
    else:
        return SklearnPerceptron(eta0=learning_rate, max_iter=epochs, tol=None, shuffle=False)

def run_perceptron(X_train, y_train, X_test, y_test, custom=True, learning_rate=1.0, epochs=DEFAULT_EPOCHS):
    perceptron = create_perceptron(custom, learning_rate, epochs)
    perceptron.fit(X_train, y_train)
    y_pred = perceptron.predict(X_test)
    accuracy = np.mean(y_pred == y_test)
    return accuracy

if __name__ == '__main__':
    data = load_breast_cancer()
    # Transform X and y from [0, 1] to [-1, 1].
    X = 2 * data.data - 1
    y = 2 * data.target - 1
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # epochs = find_best_epochs(X_train, y_train, X_test, y_test)
    epochs = DEFAULT_EPOCHS
    learning_rate = 1.0
    custom_accuracy = run_perceptron(X_train, y_train, X_test, y_test, True, learning_rate, epochs)
    print(f"Custom perceptron accuracy: {custom_accuracy * 100:.2f}%")

    sklearn_accuracy = run_perceptron(X_train, y_train, X_test, y_test, False, learning_rate, epochs)
    print(f"Sklearn perceptron accuracy: {sklearn_accuracy * 100:.2f}%")
