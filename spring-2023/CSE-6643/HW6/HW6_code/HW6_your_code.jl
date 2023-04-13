using LinearAlgebra

#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix T and modifies it 
# in place to Hessenberg form using Householder reduction.
function hessenberg_form!(T)
    m = size(T, 1)
    for k = 1:m-2
        v = T[k+1:m, k]
        v[1] += norm(v) * sign(v[1]) 
        v ./= norm(v)

        # I implemented Trefethen & Bau 26.1 without allocations, but I could not get it to preserve the eigenvalues,
        # so I implemented the approach in
        # http://pi.math.cornell.edu/~web6140/TopTenAlgorithms/QRalgorithm.html.
        # I have partially unrolled the array version presented there, but there is still the very inefficient section at the bottom.
        # Also there's the implicit array copy for `v` above. Eeps.
        # Much more work to do here, sorry I didn't get this completely done.
        # Array version:
        #     Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        #     T[k+1:m,k] = Q1*T[k+1:m,k]
        #     T[k,k+1:m] = Q1*T[k,k+1:m]
        #     T[k+1:m,k+1:m] = Q1*T[k+1:m,k+1:m]*Q1'

        # Loop version:
        for i = k+1:m
            sum = 0
            for j = k+1:m
                sum += ((i == j ? 1 : 0) - 2 * v[i-k] * v[j-k]) * T[j, k]
            end
            T[i, k] = sum
        end
        for j = k+1:m
            sum = 0
            for i = k+1:m
                sum += ((i == j ? 1 : 0) - 2 * v[i-k] * v[j-k]) * T[k, i]
            end
            T[k, j] = sum
        end

        # Just need to unroll these last part two multiplies
        Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        T[k+1:m,k+1:m] *= Q1'
        T[k+1:m,k+1:m] = Q1 * T[k+1:m,k+1:m]
    end

    # Non-tri-diagonal values are close to zero.
    # Here, we set them to exactly zero.
    for i = 1:m
        for j = i+2:m
            T[i,j] = T[j,i] = 0
        end
    end
end

#----------------------------------------
# Problem b
#----------------------------------------
# This funciton takes in a matrix T in Hessenberg form
# and runs a single iteration of the unshifted QR Algorithm 
# using Givens rotations
function givens_qr!(T)
    m = size(T, 1)

    # Store Givens rotations for use in second loop.
    G = Vector{Tuple{Float64, Float64}}(undef, m-1)
    for k = 1:m-1
        # Compute and store Givens rotation.
        a, b = T[k, k], T[k + 1, k]
        r = sqrt(a^2 + b^2)
        c, s = a / r, -b / r
        G[k] = (c, s)
        for j = 1:m
            τ₁ = T[k, j]
            τ₂ = T[k + 1, j]
            T[k, j] = τ₁ * c - τ₂ * s
            T[k + 1, j] = τ₁ * s + τ₂ * c
        end
    end

    for k = 1:m-1
        c, s = G[k]
        for j = 1:m
            τ₁ = T[j, k]
            τ₂ = T[j, k + 1]
            T[j, k] = τ₁ * c - τ₂ * s
            T[j, k + 1] = τ₁ * s + τ₂ * c
        end
    end
end

# Wilkinson shift for symmetric matrices.
function wilkinson_shift(a, b, c)
    δ = (a-c)/2
    return c - sign(δ)*b^2/(abs(δ) + sqrt(δ^2+b^2))
end

#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix T in Hessenberg form and 
# implements the practical QR algorithm with shifts. 
# The input shift dictates which shift type your 
# algorithm should use. For shift = "single" implement the single shift 
# and for shift = "wilkinson" implement the Wilkinson shift
function practical_QR_with_shifts!(T, shift)
    m = size(T, 1)
    m == 1 && return # Base case. All eigenvalues found.

    ϵ = 1e-11 # Tolerance for convergence.
    max_iter = 100
    while abs(T[m, m-1]) > ϵ && (max_iter -= 1) > 0
        if shift == "single"
            μ = T[m, m]
        elseif shift == "wilkinson"
            μ = wilkinson_shift(T[m-1,m-1], T[m,m], T[m-1,m])
        end

        for i = 1:m
            T[i, i] -= μ
        end
        givens_qr!(T)
        for i = 1:m
            T[i, i] += μ
        end
    end
    practical_QR_with_shifts!(@view(T[1:m-1, 1:m-1]), shift)
end

#----------------------------------------
# Problem d
#----------------------------------------
