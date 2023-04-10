using LinearAlgebra

#----------------------------------------
# Problem a
#----------------------------------------
# This function takes in a matrix T and modifies it 
# in place to Hessenberg form using Householder reduction.
# function hessenberg_form!(T)
#     m = size(T, 1)

#     for k = 1:m-2
#         v = @view(T[k+1:m,k])
#         v1 = sign(v[1]) * norm(v)
#         v[1] += v1
#         nv = norm(v)
#         v = copy(v) # I can't find a way to avoid this copy.
#         for j = k+1:m
#             α = v' * @view(T[k+1:m,j]) / nv
#             for i = k+1:m
#                 T[i,j] -= 2α * v[i-k] / nv
#             end
#         end

#         T[k+1,k] = -v1
#         T[k+2:m,k] .= 0

#         for j = k+1:m
#             β = v' * @view(T[k+1:m, j]) / nv
#             for i = 1:m
#                 T[i, j] -= 2β * v[j-k] / nv
#             end
#         end
#     end
# end

function hessenberg_form!(T)
    m = size(T, 1)
    for k = 1:m-2
        v = T[k+1:m, k]
        v[1] += norm(v) * sign(v[1]) 
        v ./= norm(v)

        # Array version:
        # Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        # T[k+1:m,k] = Q1*T[k+1:m,k]
        # T[k,k+1:m] = Q1*T[k,k+1:m]
        # T[k+1:m,k+1:m] = Q1*T[k+1:m,k+1:m]*Q1'

        # Loop version:
        # Equivalent to:
        # ```
        # Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        # T[k+1:m,k] = Q1*T[k+1:m,k]
        # ```
        for i = 1:m-k
            sum = 0
            for j = 1:m-k
                sum += ((i == j ? 1 : 0) - 2 * v[i] * v[j]) * T[k+j, k]
            end
            T[k+i, k] = sum
        end

        # Equivalent to:
        # ```
        # Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        # T[k,k+1:m] = Q1*T[k,k+1:m]
        # ```
        for j = 1:m-k
            sum = 0
            for i = 1:m-k
                sum += ((i == j ? 1 : 0) - 2 * v[i] * v[j]) * T[k, k+i]
            end
            T[k, k+j] = sum
        end

        Q1 = Matrix{Float64}(I, m-k, m-k) - 2(v*v')
        T[k+1:m,k+1:m] = Q1*T[k+1:m,k+1:m]*Q1'

        # Loop version, avoiding creation of identity matrix:
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

# eye = Matrix{Float64}(I, m, m)
#----------------------------------------
# Problem c
#----------------------------------------
# This function takes in a matrix T in Hessenberg form and 
# implements the practical QR algorithm with shifts. 
# The input shift dictates which shift type your 
# algorithm should use. For shift = "single" implement the single shift 
# and for shift = "wilkinson" implement the Wilkinson shift

# Wilkinson shift for symmetric matrices.
function wilkinson_shift(a, b, c)
    δ = (a-c)/2
    return c - sign(δ)*b^2/(abs(δ) + sqrt(δ^2+b^2))
end

function practical_QR_with_shifts!(T, shift)
    m = size(T, 1)
    ϵ = 1e-12
    max_iterations = 100 * m
    for _ = 1:max_iterations
        for k = 1:m-1
            # Check for deflation
            if abs(T[k+1, k]) < ϵ * (abs(T[k, k]) + abs(T[k+1, k+1]))
                T[k+1, k] = 0
                k += 1
                continue
            end

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

        # Check for convergence
        sub_diag_norm = 0
        for i = 1:m-1
            sub_diag_norm += abs(T[i+1, i])
        end

        if sub_diag_norm < ϵ
            break
        end
    end
end

# From http://pi.math.cornell.edu/~web6140/TopTenAlgorithms/QRalgorithm.html
# function QRwithShifts( A::Matrix )
#     n = size(A,1)
#     myeigs = zeros(n)
#     if ( n == 1 )
#         myeigs[1] = A[1,1]
#     else
#         I = eye( n )
#         # Reduction to Hessenberg form:
#         A = HessenbergReduction( A )
#         # Let's start the shifted QR algorithm with 
#         while( norm(A[n,n-1]) > 1e-10 )
#             mu = wilkinson_shift( A[n-1,n-1], A[n,n], A[n-1,n] )
#             # This line should use faster Hessenberg reduction:
#             (Q,R) = qr(A - mu*I)
#             # This line needs speeding up, currently O(n^3) operations!: 
#             A = R*Q + mu*I
#         end
#         # Deflation and recurse:
#         myeigs = [A[n,n] ; QRwithShifts( A[1:n-1, 1:n-1] )]
#     end
#     return myeigs
# end

#----------------------------------------
# Problem d
#----------------------------------------


