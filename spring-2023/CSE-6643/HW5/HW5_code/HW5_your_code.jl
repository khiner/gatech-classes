import Pkg.instantiate
instantiate()
using LinearAlgebra
using CairoMakie

"""
Create a diagonalizable matrix of size m x m, with the provided eigenvalues Λ. 
Performs a QR decomposition of a random matrix to create an orthogonal matrix Q, and returns QΛQ'.
"""
function random_matrix(m, Λ)
    Q, _ = qr(randn(m, m))
    return Q * Diagonal(Λ) * Q'
end

function orthogonal_iteration(A, k_max)
    # Initialize the matrix Q to the identity matrix.
    Q = Matrix{Float64}(I, size(A))

    for k = 1:k_max
        Q, R = qr(A * Q)
        println(round.(diag(R), digits=4))
    end
end

# Initialize the matrix.
m = 8
A = random_matrix(m, 3. .^(0:m-1))

orthogonal_iteration(A, 10)
