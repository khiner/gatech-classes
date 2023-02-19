#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in the LU factorization
# together with the permutation P (in the form of
# an array of integers such that P[i] = j means that multiplication with P moves the i=th row to the j-th position)
# It should modify the input variable x in place
function substitution!(x, LU, P = 1:size(x, 1))
    m = size(x, 1)
    permute!(x, P)
    # Forward substitution
    for k = 1:m-1
        for i = k+1:m
            x[i] -= LU[i,k] * x[k]
        end
    end
    # Backward substitution
    for k = m:-1:1
        x[k] /= LU[k,k]
        for i = 1:k-1
            x[i] -= LU[i,k] * x[k]
        end
    end
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and modifies
# it in place, such that is contains the stricly
# lower triangular part of L together with the
# upper triangular part of U.
function unpivoted_LU!(A)
    m = size(A, 1)
    for k = 1:m-1
        for i = k+1:m
            A[i,k] /= A[k,k]
            for j = k+1:m
                A[i,j] -= A[i,k] * A[k,j]
            end
        end
    end
end

#----------------------------------------
# Problem d
#----------------------------------------
# This function takes in a matrix A and modifies
# it in place, such that is contains the stricly
# lower triangular part of L together with the
# upper triangular part of U.
# It uses row-pivoting and stores the resulting 
# row permutation in the array P
function pivoted_LU!(A)
    m = size(A, 1)
    P = collect(1:m)
    for k = 1:m-1
        # Find pivot and update permutation
        pivot_idx = argmax(abs.(A[k:m,k]))
        pivot_idx += k - 1
        if pivot_idx != k
            A[[k,pivot_idx],:] = A[[pivot_idx,k],:]
            P[[k,pivot_idx]] = P[[pivot_idx,k]]
        end
        # Matrix factorization
        for i = k+1:m
            A[i,k] /= A[k,k]
            for j = k+1:m
                A[i,j] -= A[i,k] * A[k,j]
            end
        end
    end
    return P
end

#----------------------------------------
# Problem e
#----------------------------------------
# Creates an m × m matrix with a particularly 
# large growth factor
function growth_matrix(m)
    A = zeros(m, m)
    for i = 1:m
        for j = 1:m
            A[i,j] = i == j ? 1 : i > j ? 10^(-6) : 10^6
        end
    end
    return A
end
