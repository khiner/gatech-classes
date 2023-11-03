import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split

from perceptron import Perceptron, DEFAULT_EPOCHS

class AdaBoost:
    def __init__(self, T=50):
        self.T = T
        self.alphas = []
        self.models = []
        self.train_accuracy_history = []
        self.test_accuracy_history = []
        self.margin_history = []

    def fit(self, X, y, X_test, y_test):
        N = len(y)
        D = np.ones(N) / N
        misclassification_ft = []
        misclassification_ft_minus_1 = []
        for t in range(self.T):
            h_t = Perceptron(epochs=DEFAULT_EPOCHS*2)
            h_t.fit(X, y, sample_weight=D)

            y_pred = h_t.predict(X)
            error_t = np.sum(D * (y_pred != y))
            alpha_t = 0.5 * np.log((1 - error_t) / max(error_t, 1e-16))

            D *= np.exp(-alpha_t * y * y_pred)
            D /= D.sum()

            self.alphas.append(alpha_t)
            self.models.append(h_t)

            self.train_accuracy_history.append(np.mean(self.predict(X) == y))
            self.test_accuracy_history.append(np.mean(self.predict(X_test) == y_test))

            # Calculate the confidence margin for this iteration and store it.
            H_x = np.zeros(X.shape[0])
            for alpha, h in zip(self.alphas, self.models):
                H_x += alpha * h.predict(X)
            margin = np.min(np.abs(H_x))
            self.margin_history.append(margin)

            # 1. Generate samples from D
            indices = np.random.choice(N, size=N, p=D)
            X_sample, y_sample = X[indices], y[indices]

            # 2. Calculate the percentage of misclassified points
            if t > 0:
                y_pred_t_minus_1 = self.models[-2].predict(X_sample)
                error_t_minus_1 = np.mean(y_pred_t_minus_1 != y_sample)
                misclassification_ft_minus_1.append(error_t_minus_1)
                
            y_pred_t = self.models[-1].predict(X_sample)
            error_t = np.mean(y_pred_t != y_sample)
            misclassification_ft.append(error_t)
        
        return misclassification_ft, misclassification_ft_minus_1

    def predict(self, X_in):
        H_x = np.zeros(X_in.shape[0])
        for alpha, h in zip(self.alphas, self.models):
            H_x += alpha * h.predict(X_in)
        return np.sign(H_x)

if __name__ == '__main__':
    data = load_breast_cancer()
    # Transform X and y from [0, 1] to [-1, 1].
    X_g = 2 * data.data - 1
    y_g = 2 * data.target - 1
    X_train_g, X_test_g, y_train_g, y_test_g = train_test_split(X_g, y_g, test_size=0.2, random_state=42)

    ada = AdaBoost(T=30)
    misclass_ft, misclass_ft_minus_1 = ada.fit(X_train_g, y_train_g, X_test_g, y_test_g)
    y_pred = ada.predict(X_test_g)
    accuracy = np.mean(y_pred == y_test_g)
    print(f"AdaBoost accuracy: {accuracy * 100:.2f}%")

    # Plotting the confidence margin
    # plt.plot(range(1, len(ada.margin_history) + 1), ada.margin_history)
    # plt.xlabel('Iterations')
    # plt.ylabel('Confidence Margin')
    # plt.title('Confidence Margin vs Number of Iterations')
    # plt.savefig('adaboost_margin.png')

    # plt.figure()
    # print(f"ta: {ada.train_accuracy_history}")
    # plt.plot(ada.train_accuracy_history, label='Training Accuracy')
    # plt.plot(ada.test_accuracy_history, label='Test Accuracy')
    # plt.xlabel('Iterations')
    # plt.ylabel('Accuracy')
    # plt.title('Accuracy vs Number of Iterations')
    # plt.legend()
    # plt.savefig('adaboost_train_test_accuracy.png')

    plt.figure()
    plt.plot(misclass_ft, label='$f_t$ misclassification rate')
    plt.plot(misclass_ft_minus_1, label='$f_{t-1}$ misclassification rate')
    plt.xlabel('Iterations')
    plt.ylabel('Misclassification Rate')
    plt.title('Misclassification Rate vs Number of Iterations')
    plt.legend()
    plt.savefig('adaboost_misclassification_rate.png')