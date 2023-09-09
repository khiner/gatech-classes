import scipy.io
import numpy as np
from myRecommender import run

cell = scipy.io.loadmat('movie_data.mat')
rate_mat = cell['train']
test_mat = cell['test']

def rmse(u, v, mat):
    mask = mat > 0
    res = np.sum(((u.dot(v.T) - mat) * mask) ** 2) / float(np.sum(mask))
    return np.sqrt(res)

lr = 5 # "Low-rank". Using the highest value from HW1.
max_iter = 1000 # For each run.

best_learning_rate = 0
best_reg_coef = 0
best_rmse = float('inf')
for learning_rate in [0.0001, 0.0002, 0.0004, 0.0005]:
    for reg_coef in [0.001, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4]:
        U, V = run(rate_mat, lr, True, learning_rate, reg_coef, max_iter)
        rmse_val = rmse(U, V, test_mat)
        print('Learning rate: {}, Reg coef: {}, RMSE: {}'.format(learning_rate, reg_coef, rmse_val))
        if rmse_val < best_rmse:
            best_rmse = rmse_val
            best_learning_rate = learning_rate
            best_reg_coef = reg_coef

print('Learning rate: {}, Reg coef: {}'.format(best_learning_rate, best_reg_coef)) 
