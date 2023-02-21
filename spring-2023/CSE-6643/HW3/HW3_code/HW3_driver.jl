import Pkg.instantiate
instantiate()
using BenchmarkTools: @ballocated
using Random: randperm
using LinearAlgebra: I, norm
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
@assert allocated_memory < 450

#----------------------------------------
# Problem c
#----------------------------------------
function plot_error_and_growth(pivoted = false, use_growth_matrix = false)
    # Returns the maximum absolute value of all elements
    # in the upper-triangular part of the matrix `A`
    # (including the diagonal).
    function max_abs_u(A)
        m = size(A, 1)
        max_abs = 0.0 # Running max value
        for i = 1:m
            for j = i:m
                max_abs = max(max_abs, abs(A[i, j]))
            end
        end
        return max_abs
    end

    ms = map(p -> 2^p, 4:11)
    growth_factors = Vector{Float64}()
    relative_errors = Vector{Float64}()
    for m in ms
        A = use_growth_matrix ? growth_matrix(m) : (randn(m, m) + I)
        b = randn(m)
        reference_x = A \ b
        A_max = maximum(abs.(A))
        if (pivoted)
            P = pivoted_LU!(A)
            substitution!(b, A, P)
        else
            unpivoted_LU!(A) 
            substitution!(b, A)
        end
        append!(growth_factors, max_abs_u(A) / A_max)
        append!(relative_errors, norm(b - reference_x) / norm(reference_x))
    end

    fig = Figure()
    
    relative_error_ax = Axis(fig[1, 1],
        title = L"Relative error for matrix with size $m$",
        xlabel = L"Matrix size $m$",
        ylabel = "Relative error",
        xscale = log2,
        xticks = ms,
    )
    scatter!(relative_error_ax, ms, relative_errors)
    # Only plotting error for growth matrix
    if !use_growth_matrix
        growth_ax = Axis(fig[2, 1],
        title = L"Growth factor for matrix with size $m$",
        xlabel = L"Matrix size $m$",
        ylabel = L"Growth factor $\rho$",
            xscale = log2,
            xticks = ms,
        )
        scatter!(growth_ax, ms, growth_factors)
    end
    save("rel_error_and_growth_" *
        (pivoted ? "pivot" : "no_pivot") *
        (use_growth_matrix ? "_gm" : "") *
        ".png",
        fig
    )
end

plot_error_and_growth(false)

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

plot_error_and_growth(true)

#----------------------------------------
# Problem e
#----------------------------------------
plot_error_and_growth(false, true)
