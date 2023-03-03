#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function classical_gram_schmidt(A)
    m, n = size(A)
    Q = zeros(m, n)
    R = zeros(n, n)

    for j = 1:n
        v = A[:,j]
        for i = 1:j-1
            R[i,j] = Q[:,i]' * A[:,j]
            v -= R[i,j] * Q[:,i]
        end
        R[j,j] = norm(v)
        Q[:,j] = v / R[j,j]
    end
    return Q, R
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function modified_gram_schmidt(A)
    m, n = size(A)
    Q = zeros(m, n)
    R = zeros(n, n)
    V = copy(A)

    for j = 1:n
        R[1:j-1,j] = Q[:,1:j-1]'*A[:,j]
        V[:,j] -= Q[:,1:j-1] * R[1:j-1,j]
        R[j,j] = norm(V[:,j])
        Q[:,j] = V[:,j] / R[j,j]
    end
    return Q, R
end

#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix A 
# and computes its QR factorization in place,
# using householder reflections.
# # It should not allocate any memory.
function householder_QR!(A)
    m, n = size(A)
    for k = 1:n
        x_norm = norm(view(A, k:m, k))
        v = A[k:m, k]
        v[1] += sign(A[k,k]) * x_norm
        v /= v[1]
        v_dot = v' * v

        # Perform `A[k:m, k:n] -= 2 * v * (v' * A[k:m, k:n])`,
        # one column at a time to avoid slicing, which allocates memory.
        for j = k:n
            inner_product = v' * view(A, k:m, j) / v_dot
            for i = k:m
                A[i,j] -= 2 * v[i-k+1] * inner_product
            end
        end
        # Perform `A[k+1:m, k] = v[2:end]`, without slicing.
        for i = k+1:m
            A[i,k] = v[i-k+1]
        end
    end
end

#----------------------------------------
# Problem d
#----------------------------------------
# These two functions take in the housholder
# QR factorization from part c and multiply them
# to a vector (mul) or solve the least squares 
# problem in A (div), in place.
# They should not allocate any memory and instead
# use the preallocated output vector to record the result. 
function householder_QR_mul!(out, x, QR)
    m, n = size(QR)
    out .= x
    for j = 1:n
        v = zeros(m)
        v[j] = 1
        v[j+1:m] = QR[j+1:m, j]
        out[j:m] -= 2 * v[j:m] * (v[j:m]' * out[j:m])
    end
end

function householder_QR_div!(out, b, QR)
    m, n = size(QR)
    out .= b
    for j = n:-1:1
        v = zeros(m)
        v[j] = 1
        v[j+1:m] = QR[j+1:m, j]
        out[j] -= 2 * v[j:m] * (v[j:m]' * out[j:m])
        out[j] /= QR[j,j]
    end
end

# Problem c
# function householder_QR!(A)
#     m, n = size(A)
#     for k = 1:min(m-1, n)
#         x = A[k:m,k]
#         e = zeros(length(x))
#         e[1] = 1
#         v = sign(x[1]) * norm(x) * e + x
#         v /= norm(v)
#         A[k:m,k:n] -= 2*v*(v'*A[k:m,k:n])
#     end
# end
