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
end

# Implements a recursive, cache oblivious algorithm
# complete this skeleton
function oblivious_add_to_A_B_times_C!(A, B, C, bks)
    i_size = size(A, 1)
    j_size = size(C, 2)
    k_size = size(B, 2)

    # If we want to further subdivide
    if min(i_size, j_size, k_size) > bks
    #If we are ready to break the recursion
    else
        add_to_A_B_times_C!(A, B, C)
    end
end

