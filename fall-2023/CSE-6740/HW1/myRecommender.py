import numpy as np

def rmse(U, V, M):
    mask = M > 0
    res = np.sum(((U.dot(V.T) - M) * mask) ** 2) / float(np.sum(mask))
    return np.sqrt(res)

def grad_U(U, V, M, reg_coef, with_reg):
    mask = M > 0
    grad = -2 * np.dot((M - np.dot(U, V.T)) * mask, V) + (2 * reg_coef * U if with_reg else 0)
    return grad

def grad_V(U, V, M, reg_coef, with_reg):
    mask = M > 0
    grad = -2 * np.dot(((M - np.dot(U, V.T)) * mask).T, U) + (2 * reg_coef * V if with_reg else 0)
    return grad

def my_recommender(rate_mat, lr, with_reg):
    """

    :param rate_mat:
    :param lr:
    :param with_reg:
        boolean flag, set true for using regularization and false otherwise
    :return:
    """

    # TODO adjust using cross validation
    max_iter = 1000
    learning_rate = 0.0004
    reg_coef = 0.01
    n_user, n_item = rate_mat.shape[0], rate_mat.shape[1]

    U = np.random.rand(n_user, lr) / lr
    V = np.random.rand(n_item, lr) / lr

    for i in range(max_iter):
        U -= learning_rate * grad_U(U, V, rate_mat, learning_rate, with_reg)
        V -= learning_rate * grad_V(U, V, rate_mat, reg_coef, with_reg)

    return U, V
