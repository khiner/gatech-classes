#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in the LU factorization
# together with the permutation P (in the form of
# an array of integers such that P[i] = j means that multiplication with P moves the i=th row to the j-th position)
# It should modify the input variable x in place
function substitution!(x, LU, P = 1 : size(x, 1)) 
   # YOUR CODE HERE
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and modifies
# it in place, such that is contains the stricly
# lower triangular part of L together with the
# upper triangular part of U.
function unpivoted_LU!(A)
    # YOUR CODE HERE
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
    # The array that will be used to keep track of the permutation
    P = collect(1 : size(A, 1))
    # YOUR CODE HERE
    # returns the array representing the permutation
    return P
end

#----------------------------------------
# Problem e
#----------------------------------------
# Creates an m × m matrix with a particularly 
# large growth factor
function growth_matrix(m)
    # YOUR CODE GOES HERE
end
