# This macro helps optimize the innermost loop in the kernel on your machine. You should not need to use it anywhere else
using LoopVectorization: @turbo
# This function takes in matrices ABC and adds B times C to the matrix A
function add_to_A_B_times_C!(A, B, C)
    @turbo for j in axes(C, 2)
        for k in axes(B, 2)
            for i in axes(A, 1)
                A[i, j] += B[i, k] * C[k, j]
            end
        end
    end
end

# This function takes in matrices ABC and adds B times C to the matrix A
# It uses blocking into blocks of size bks
# Make sure that your function does not allocate memory
function add_to_A_B_times_C!(A, B, C, bks)
    for i in range(1, size(A, 1) ÷ bks + 1)
        for j in range(1, size(A, 2) ÷ bks + 1)
            i_rng = (i - 1)*bks + 1:min(i*bks, size(A, 1))
            j_rng = (j - 1)*bks + 1:min(j*bks, size(A, 2))
            add_to_A_B_times_C!(@view(A[i_rng, j_rng]), @view(B[i_rng, :]), @view(C[:, j_rng]))
        end
    end
end

# Implements a recursive, cache oblivious algorithm
# complete this skeleton
# Obtain an algorithm that equally divides each of its input matrices into four parts
# and then recurses on the resulting eight subproblems,
# until one of the problems reaches a size below bks.
function oblivious_add_to_A_B_times_C!(A, B, C, bks)
    i_size = size(A, 1)
    j_size = size(C, 2)
    k_size = size(B, 2)

    # If we want to further subdivide
    if min(i_size, j_size, k_size) > bks
        for i in range(1,2)
            for j in range(1,2)
                i_rng = (i - 1)*(i_size ÷ 2) + 1:(i == 1 ? i_size ÷ 2 : i_size)
                j_rng = (j - 1)*(j_size ÷ 2) + 1:(j == 1 ? j_size ÷ 2 : j_size)
                oblivious_add_to_A_B_times_C!(@view(A[i_rng, j_rng]), @view(B[i_rng, :]), @view(C[:, j_rng]), bks)
            end
        end
    #If we are ready to break the recursion
    else
        add_to_A_B_times_C!(A, B, C)
    end
end

