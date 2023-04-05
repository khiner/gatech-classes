using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, triu, tril, tr, diagm, diag, qr
using CairoMakie
include("HW6_your_code.jl")


#----------------------------------------
# Problem a 
#----------------------------------------
########################################
m = 20
T = rand(m,m)
T = T'T
traceA = tr(T) 
hessenberg_form!(T)
@assert sum(tril(T,-2)) ≈ 0
@assert tr(T) ≈ traceA
T = randn(5, 5)
#allocated_memory = @ballocated  hessenberg_form!(T)
#@assert allocated_memory == 0
println("Passed part (a) test")


#----------------------------------------
# Problem b 
#----------------------------------------
########################################
m = 20
T = rand(m,m)
T = T'T
hessenberg_form!(T)
Q,R =qr(T)
RQ = R*Q
givens_qr!(T)
@assert abs.(T) ≈ abs.(RQ)
println("Passed part (b) test")

#----------------------------------------
# Problem c
#----------------------------------------
########################################
m = 10
A = rand(m,m)
Q,_ = qr(A)
λ = 3 .^ range(0,m-1)
Σ = diagm(λ)
T = Q'Σ*Q
hessenberg_form!(T)
shift = "wilkinson"
practical_QR_with_shifts!(T,shift)
@assert λ ≈ sort(diag(T))


#----------------------------------------
# Problem d
#----------------------------------------
# YOUR CODE GOES HERE

