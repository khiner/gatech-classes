#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function classical_gram_schmidt(A)
    # YOUR CODE HERE
    return Q, R
end

#----------------------------------------
# Problem b
#----------------------------------------
# This function takes in a matrix A and returns 
# a reduced QR factorization with factors Q and R.
# It should not modify A
function modified_gram_schmidt(A)
    # YOUR CODE HERE
    return Q, R
end

#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix A 
# and computes its QR factorization in place,
# using householder reflections.
# It should not allocate any memory.  
function householder_QR!(A)
    # YOUR CODE HERE
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
    # YOUR CODE HERE
end

function householder_QR_div!(out, b, QR)
    # YOUR CODE HERE
end