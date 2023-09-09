import numpy as np

def rmse(u, v, mat):
    mask = mat > 0
    res = np.sum(((u.dot(v.T) - mat) * mask) ** 2) / float(np.sum(mask))
    return np.sqrt(res)

def grad_U(U, V, M, reg_coef, with_reg):
    mask = M > 0
    return -2 * np.dot((M - np.dot(U, V.T)) * mask, V) + (2 * reg_coef * U if with_reg else 0)

def grad_V(U, V, M, reg_coef, with_reg):
    mask = M > 0
    return -2 * np.dot(((M - np.dot(U, V.T)) * mask).T, U) + (2 * reg_coef * V if with_reg else 0)

def run(M, lr, with_reg = True, learning_rate = 0.0001, reg_coef = 0.01, max_iter = 1000):
    n_user, n_item = M.shape[0], M.shape[1]
    U = np.random.rand(n_user, lr) / lr
    V = np.random.rand(n_item, lr) / lr

    error = rmse(U, V, M)
    for _ in range(max_iter):
        U -= learning_rate * grad_U(U, V, M, reg_coef, with_reg)
        V -= learning_rate * grad_V(U, V, M, reg_coef, with_reg)
        new_error = rmse(U, V, M)
        # Error delta found empirically to be similar to stopping when validation error increases.
        if abs(new_error - error) < 1e-4:
            break
        error = new_error

    return U, V

def my_recommender(rate_mat, lr, with_reg):
    """

    :param rate_mat:
    :param lr:
    :param with_reg:
        boolean flag, set true for using regularization and false otherwise
    :return:
    """

    max_iter = 1500
    learning_rate = 0.0001
    reg_coef = 0.2
    return run(rate_mat, lr, with_reg, learning_rate, reg_coef, max_iter)
