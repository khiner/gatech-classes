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
allocated_memory = @ballocated hessenberg_form!(T)
# @assert allocated_memory = 0
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
println("Passed part (c) test")

#----------------------------------------
# Problem d
#----------------------------------------
# YOUR CODE GOES HERE

# You just need to show the convergence for one of the eigenvalues. It is up to you which one you select. 
# Use convergence criteria similar to HW5 P3.

"""
Plot the convergence of the estimated eigenvalues to the true eigenvalues.
  * The first subplot shows the estimated and true eigenvalues.
  * The second subplot shows the actual and theoretical convergence rates.
Arguments:
  * `λ_history` is a dictionary mapping shift type to a vector of the estimated eigenvalues at each iteration.
  * `Λ_true` are the true eigenvalues.
  * `outfile` is the path of the file to save the plot to.
"""
function plot_convergence(λ_history::Dict{String, Vector{Float64}}, Λ_true, λ_index, outfile="convergence.png")
    fig = Figure()
    shifts = λ_history |> keys |> collect |> sort
    values_ax = Axis(fig[1, 1],
        title = L"Estimated vs. true $\lambda$",
        xlabel = L"Iteration $k$",
        ylabel = "Estimated eigenvalue",
        yscale = Makie.Symlog10(3.0), # Symlog scale allows for negative values.
        yticks = [-3 .^ (8:-1:1); 0; 3 .^ (1:8)],
        ytickformat = "{:d}",
    )
    conv_ax = Axis(fig[2, 1],
        title = L"Actual vs. theoretical convergence of $\lambda$",
        xlabel = L"Iteration $k$",
        ylabel = "Convergence",
        yscale = log10,
    )

    linestyles = [:solid, :dash]
    colors = map(i ->  Makie.wong_colors()[i], 1:length(shifts))
    for (i, (shift, history)) in enumerate(λ_history)
        n = length(history)
        color = colors[i]

        # Plot the estimated and true eigenvalues.
        lines!(values_ax, 1:n, history, linestyle=linestyles[1], color=color)
        hlines!(values_ax, [Λ_true[λ_index]], linestyle=linestyles[2], color=color)
 
        # Plot the actual and theoretical convergence rates.
        conv = abs.(history .- Λ_true[λ_index]) ./ Λ_true[λ_index]
        theoretical_conv_base =
            λ_index == 1 ? abs(Λ_true[2] / Λ_true[1]) :
            λ_index == m ? abs(Λ_true[m] / Λ_true[m - 1]) :
            max(abs(Λ_true[λ_index+1] / Λ_true[λ_index]), abs(Λ_true[λ_index] / Λ_true[λ_index-1]))
        theoretical_conv = theoretical_conv_base .^ (1:n)
        lines!(conv_ax, 1:n, conv, linestyle=linestyles[1], color=color)
        lines!(conv_ax, 1:n, theoretical_conv, linestyle=linestyles[2], color=color)
    end

    group_style = [
        LineElement(linestyle = ls, color = :black)
        for ls in linestyles
    ]
    group_color = [
        PolyElement(color = color, strokecolor = :transparent)
        for color in colors
    ]
    Legend(fig[1, 2],
        [group_style, group_color],
        [
            ["Estimated", "True"], shifts,
        ],
        ["Shift", "Eigenvalue"]
    )
    Legend(fig[2, 2],
        [group_style, group_color],
        [
            ["Actual", "Theoretical"], shifts,
        ],
        ["Shift", "Eigenvalue"]
    )

    save(outfile, fig)
end


"""
Same algorithm as `practical_QR_with_shifts!`, but plots the convergence
of the eigenvalues for different shifts.
`Λ_true` are the true eigenvalues.
Plots the convergence of the eigenvalue value at `λ_index`.
"""
function monitored_practical_QR_with_shifts(T, Λ_true, λ_index)
    shifts = ["single", "wilkinson"]
    λ_history = Dict(shift => Float64[] for shift in shifts)
    T_orig = copy(T)

    m = size(T, 1)
    # m == 1 && return # Base case. All eigenvalues found.

    ϵ = 1e-11 # Tolerance for convergence.
    max_iter = 100
    for shift in shifts
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
 
            # Monitoring:
            λ_est = T[λ_index, λ_index]
            push!(λ_history[shift], λ_est)
        end
        T .= T_orig
    end
    # practical_QR_with_shifts!(@view(T[1:m-1, 1:m-1]), shift)

    plot_convergence(λ_history, Λ_true, λ_index)
end

# Random.seed!(1234)
m = 10
A = rand(m,m)
Q,_ = qr(A)
Λ_true = 3. .^ range(0,m-1)
Σ = diagm(Λ_true)
A = Q'Σ*Q
hessenberg_form!(A)
λ_index = m
monitored_practical_QR_with_shifts(A, Λ_true, λ_index)
