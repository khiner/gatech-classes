import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, istriu, triu, qr, Diagonal
using CairoMakie
include("HW4_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
A = randn(30, 20) 
b = randn(30)
Q, R = classical_gram_schmidt(A) 
@assert Q' * Q ≈ I
@assert Q * R ≈ A

#----------------------------------------
# Problem b 
#----------------------------------------
########################################
A = randn(30, 20) 
b = randn(30)
Q, R = modified_gram_schmidt(A) 
@assert Q' * Q ≈ I
@assert Q * R ≈ A

#----------------------------------------
# Problem c
#----------------------------------------
########################################
A = randn(25, 20) 
allocated_memory = @ballocated  householder_QR!(A)
@assert allocated_memory == 0
A = randn(25, 20)
true_R = Matrix(qr(A).R)
householder_QR!(A)
# Checks if the R part of the factorization is correct
@assert vcat(true_R, zeros(5,20)) ≈ triu(A)

#----------------------------------------
# Problem d
#----------------------------------------
########################################
# Testing for memory allocation:
A = randn(25, 20) 
householder_QR!(A)
QR = A
x = randn(20)
b = randn(25)
out_mul = randn(25)
out_div = randn(20)

allocated_memory_mul = @ballocated  householder_QR_mul!(out_mul, x, QR)
allocated_memory_div = @ballocated  householder_QR_div!(out_div, b, QR)
@assert allocated_memory_mul == 0
@assert allocated_memory_div == 0

# Testing for correctness:
A = randn(25, 20) 
x = randn(20)
b = randn(25)
out_mul = randn(25)
out_div = randn(20)
true_mul = A * x 
true_div = A \ b 

householder_QR!(A)
QR = A
householder_QR_mul!(out_mul, x, QR)
householder_QR_div!(out_div, b, QR)

# checks whether the results are approximately correct
@assert true_mul ≈ out_mul
@assert true_div ≈ out_div


#----------------------------------------
# Problem e
#----------------------------------------
# YOUR CODE GOES HERE

function compare_errors()
    classical_errors = Vector{Float64}()
    modified_errors = Vector{Float64}()

    classical_r_errors = Vector{Float64}()
    modified_r_errors = Vector{Float64}()
    householder_r_errors = Vector{Float64}()
    ms = map(p -> 2^p, 4:11)
    for m = ms
        U, X = qr(randn(m)) # Set `U` to be a random orthogonal matrix.
        V, Y = qr(randn(m)) # Set `V` to be a random orthogonal matrix.
        S = Diagonal(2. .^ ((-1:(m-2)) ./ 4)) # Set `S` to be a diagonal matrix with exponentially graded entries.
        A = U*S*V # Set `A` to a matrix with these elements as singular values.

        Qc, Rc = classical_gram_schmidt(A) 
        Qm, Rm = modified_gram_schmidt(A)

        append!(classical_errors, norm(Qc * Rc - A) / norm(A))
        append!(modified_errors, norm(Qm * Rm - A) / norm(A))
        true_R = qr(A).R
        append!(classical_r_errors, norm(Rc - true_R) / norm(true_R))
        append!(modified_r_errors, norm(Rm - true_R) / norm(true_R))
        mod_A = copy(A)
        householder_QR!(mod_A)
        append!(householder_r_errors, norm(triu(mod_A) - true_R) / norm(true_R))
    end

    println(classical_errors)
    println(modified_errors)
    println(modified_r_errors)
    fig = Figure()
    recon_ax = Axis(fig[1, 1],
        title = L"Reconstruction error for matrix with size $m$",
        xlabel = L"Matrix size $m$",
        ylabel = "Relative reconstruction error",
        xscale = log2,
        xticks = ms,
        yscale = log10,
    )
    scatter!(recon_ax, ms, classical_errors, label="Classical GS")
    scatter!(recon_ax, ms, modified_errors, label="Modified GS")
    axislegend(recon_ax, position = :lt)

    r_ax = Axis(fig[2, 1],
        title = L"$R$ error for matrix with size $m$",
        xlabel = L"Matrix size $m$",
        ylabel = L"$R$ error",
        xscale = log2,
        xticks = ms,
        yscale = log10,
    )
    scatter!(r_ax, ms, classical_r_errors, label="Classical GS")
    scatter!(r_ax, ms, modified_r_errors, label="Modified GS")
    scatter!(r_ax, ms, householder_r_errors, label="Householder QR")
    axislegend(r_ax, position = :lt)

    save("error_compare.png", fig)
end

compare_errors()
