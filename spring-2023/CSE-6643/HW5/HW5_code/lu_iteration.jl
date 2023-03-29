import Pkg.instantiate
instantiate()

using LinearAlgebra

function simultaneous_iteration(A, k_max)
    T = A
    Q = Matrix{Float64}(I, size(A))
    for _ = 1:k_max
        Z = A * Q
        Q, R = qr(Z)
        T = Q' * A * Q
    end

    return T
end

function qr_iteration(A, k_max)
    T = A
    for _ = 1:k_max
        Q, R = qr(T)
        T = R * Q
    end

    return T
end

function simultaneous_lu_iteration(A, k_max)
    G = Matrix{Float64}(I, size(A))
    T = A
    for _ = 1:k_max
        G_prev = G
        Z = A * G
        G, R = lu(Z, Val(false))
        T = R * (G_prev \ G) # R * inv(G_prev) * G
        # Note: Could just use G^-1 * A * G as in simultaneous QR iteration,
        # but this is a bit more efficient, since it only involves triangular matrices.
    end

    return T
end

function lu_iteration(A, k_max)
    T = A
    for _ = 1:k_max
        L, U = lu(T, Val(false))
        T = U * L
    end

    return T
end

# `A` must be symmetric positive definite.
# Two steps of LU iteration are equivalent to one step of QR iteration.
function lr_iteration(A, k_max)
    T = A
    for _ = 1:k_max
        L = cholesky(T).L
        T = L' * L
    end

    return T
end

function rand_matrix(m, Λ)
    Q, _ = qr(randn(m, m))
    return Q * Diagonal(Λ) * Q'
end

function rand_spd_matrix(m, Λ)
    A = rand_matrix(m, sqrt.(Λ))
    return A * A'
end

m = 8
Λ_true = 3. .^(0:m-1)
A = rand_matrix(m, Λ_true)
# A = rand_spd_matrix(m, Λ_true)
# A = randn(m, m)

# println("True eigenvalues:")
# display(reverse(Λ_true))

println("Simultaneous iteration:")
display(diag(simultaneous_iteration(A, 12)))

println("QR iteration:")
display(diag(qr_iteration(A, 12)))

println("Simultaneous LU iteration:")
display(diag(simultaneous_lu_iteration(A, 24)))

println("LU iteration:")
display(diag(lu_iteration(A, 24)))

println("Est. eigenvalues:")
display(eigvals(A))

# println("LR iteration:")
# display(diag(lr_iteration(A, 24)))
