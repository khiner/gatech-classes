using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, triu, tril, tr, diagm, diag, qr, eigvals, istriu
using CairoMakie
include("HW7_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
m = 5
A = rand(m,m)
A = A * A' 
H = arnoldi(A, randn(m), m)
@assert istriu(H[2 : end, 1 : (end - 1)])
@assert eigvals(H) ≈ eigvals(A) 

println("Passed part (a) test")


#----------------------------------------
# Problem b 
#----------------------------------------
########################################
m = 5
A = rand(m,m)
A = A * A' 
α, β = lanczos(A, randn(m), m)
H = diagm(-1 => β, 0 => α, 1 => β)
@assert eigvals(H) ≈ eigvals(A)
println("Passed part (b) test")

#----------------------------------------
# Problem c
#----------------------------------------
########################################

# Your code here

