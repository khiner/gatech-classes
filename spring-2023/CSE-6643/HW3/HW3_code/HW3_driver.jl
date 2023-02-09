import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using Random: randperm
using LinearAlgebra: I
using CairoMakie
include("HW3_your_code.jl")


#----------------------------------------
# Problem a + b
#----------------------------------------
########################################
A = randn(20, 20) + 100 * I
b = randn(20)
reference_x = A \ b
unpivoted_LU!(A) 
substitution!(b, A)
@assert reference_x ≈ b

allocated_memory = @ballocated  unpivoted_LU!(A)
allocated_memory += @ballocated  substitution!(b, A)
@assert allocated_memory == 0

#----------------------------------------
# Problem c
#----------------------------------------
# YOUR CODE GOES HERE


#----------------------------------------
# Problem d
#----------------------------------------
########################################
A = (randn(20, 20) + 100 * I)[randperm(20), :]
b = randn(20)
reference_x = A \ b
P = pivoted_LU!(A) 
substitution!(b, A, P)
@assert reference_x ≈ b

#----------------------------------------
# Problem e
#----------------------------------------
# YOUR CODE GOES HERE
