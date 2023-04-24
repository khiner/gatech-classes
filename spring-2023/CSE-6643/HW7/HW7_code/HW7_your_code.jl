

#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix A and returns the 
# kmax × kmax supmatrix of itsupper Hessenberg form 
# The matrix A should be only accessed through 
# (kmax - 1) matrix-vector products
function arnoldi(A, q1, kmax)
    # Your code here
    m = kmax + 1
    H = zeros(m, m - 1)
    Q = zeros(m - 1, m)
    Q[:,1] = q1/norm(q1)
    for i = 1:kmax
        v = A*Q[:,i]
        for j = 1:i
            H[j,i] = Q[:,j]'*v
            v -= H[j,i]*Q[:,j]
        end
        H[i+1,i] = norm(v)
        Q[:,i+1] = v/H[i+1,i]
    end
    H = copy(H[1:m-1,:])
    return H
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and returns the 
# kmax × kmax supmatrix of its tridiagonal form 
# computed by the Lanczos iteration.
# The matrix A should be only accessed through 
# (kmax - 1) matrix-vector products
# The output vectors should be the diagonal (α)
# and the offdiagonal (β) of the tridiagonal matrix
function lanczos(A, q1, kmax)
    n = size(A, 1)
    α = zeros(kmax)
    β = zeros(kmax)

    Q = zeros(n, kmax + 1)
    Q[:, 1] = q1 / norm(q1)

    for k = 1:kmax
        r = A * Q[:, k]
        if k > 1
            r -= β[k - 1] * Q[:, k - 1]
        end

        α[k] = r' * Q[:, k]
        if k < kmax
            r -= α[k] * Q[:, k]
            β[k] = norm(r)
            Q[:, k + 1] = r / β[k]
        end
    end

    return α, β[1:end-1]
end
