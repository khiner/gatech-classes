import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using LinearAlgebra: I, norm, istriu, triu, qr
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
