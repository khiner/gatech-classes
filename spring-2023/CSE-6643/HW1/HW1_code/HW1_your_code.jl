# This function takes in a matrix A and a vector v and writes their product into the vector u
function u_is_A_times_v!(u, A, v)
    for i in 1:size(A, 1) # ith row of A
        u[i] = 0.0
        for j in 1:size(A, 2) # jth column of A
            u[i] += A[i,j] * v[j]
        end
    end
end

# This function takes in matrices ABC and writes B times C into the matrix A
function A_is_B_times_C!(A, B, C)
    for B_i in 1:size(B, 1) # ith row of B
        for C_j in 1:size(C, 2) # jth column of C
            A[B_i,C_j] = 0.0
            for B_j in 1:size(B, 2) # jth column of B
                A[B_i,C_j] += B[B_i,B_j] * C[B_j,C_j]
            end
        end
    end
end
