import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import sys
from myRecommender import run

cell = scipy.io.loadmat('movie_data.mat')
rate_mat = cell['train']
test_mat = cell['test']

def rmse(u, v, mat):
    mask = mat > 0
    res = np.sum(((u.dot(v.T) - mat) * mask) ** 2) / float(np.sum(mask))
    return np.sqrt(res)

max_iter = 1500
learning_rate = 0.0001
reg_coef = 0.2
lrs = [1, 3, 5, 7, 9, 11, 13, 15]
rmses = np.zeros((len(lrs), 2)) # 2 for train and test.
for i, lr in enumerate(lrs):
    U, V = run(rate_mat, lr, True, learning_rate, reg_coef, max_iter)
    rmse_val = rmse(U, V, rate_mat)
    rmse_test = rmse(U, V, test_mat)
    rmses[i] = [rmse_val, rmse_test]

plt.plot(lrs, rmses[:, 0], label='Train')
plt.plot(lrs, rmses[:, 1], label='Test')
plt.xlabel('Low-rank')
plt.ylabel('RMSE')
plt.legend()
plt.grid()
plt.title('Model performance with varying low-rank')
plt.savefig(sys.argv[1] if len(sys.argv) > 1 else 'evaluate.png')
