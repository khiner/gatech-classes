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
#  It should not allocate any memory.
function householder_QR!(A)
    m, n = size(A)
    for k = 1:n
        v = A[k:m, k]
        v[1] += sign(A[k,k]) * norm(@view(A[k:m, k]))
        v /= v[1]
        @view(A[k:m, k:n]) .-= 2 * v * (v' * @view(A[k:m, k:n])) / (v' * v)
        @view(A[k+1:m, k]) .= @view(v[2:end])
    end
end

#----------------------------------------
# Problem d
#----------------------------------------
# These two functions take in the householder
# QR factorization from part c and multiply them
# to a vector (mul) or solve the least squares 
# problem in A (div), in place.
# They should not allocate any memory and instead
# use the preallocated output vector to record the result. 

# Replace `out` with $\mtx{Q}\mtx{R}\vct{x}$.
function householder_QR_mul!(out, x, QR)
    m, n = size(QR)
 
    # Compute Rx.
    # Multiply the upper triangular part of QR by x and put it in `out`.
    for i = 1:m
        out[i] = (@view(QR[i,i:n]))' * @view(x[i:n])
    end

    # Compute QRx using algorithm from slide 7.
    for k = n:-1:1
        v = @view(QR[k:m, k])
        v1_orig = v[1]
        v[1] = 1.0
        α = 2 * v' * @view(out[k:m]) / (v' * v)
        for i = k:m
            out[i] -= v[i-k+1] * α
        end
        v[1] = v1_orig
    end
end

# Replace `out` with $\mtx{R}^{-1}\mtx{Q}^{*}\vct{b}$.
function householder_QR_div!(out, b, QR)
    m, n = size(QR)

    # Compute Q'b.
    for k = 1:n
        v = @view(QR[k:m, k])
        v1_orig = v[1]
        v[1] = 1.0
        α = 2 * (v' * @view(b[k:m])) / (v' * v)
        for i = k:m
            b[i] -= v[i-k+1] * α
        end
        v[1] = v1_orig
    end

    # Solve Rx = Q'b using backward substitution.
    for i = n:-1:1
        out[i] = (b[i] - (@view(QR[i,i+1:end]))' * @view(out[i+1:end])) / QR[i,i]
    end
end
