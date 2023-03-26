import Pkg.instantiate
instantiate()

using LinearAlgebra
using CairoMakie

"""
Create a diagonalizable matrix of size m x m, with the provided eigenvalues Λ. 
Performs a QR decomposition of a random matrix to create an orthogonal matrix Q, and returns QΛQ'.
"""
function random_matrix(m, Λ)
    Q, _ = qr(randn(m, m))
    return Q * Diagonal(Λ) * Q'
end

function plot_convergence(history::Dict{Int64, Vector{Float64}}, Λ_true::Vector{Float64}, outfile="convergence.png")
    fig = Figure()
    I = keys(history)
    I_str = join(I, ", ")
    ax = Axis(fig[1, 1],
        title = L"Convergence of $\Lambda_i$ for $i \in [%$I_str]",
        xlabel = L"Iteration $k$",
        ylabel = "Estimated eigenvalue",
        yscale = Makie.Symlog10(3.0), # Symlog scale allows for negative values.
        yticks = [-3 .^ (8:-1:0); 0; 3 .^ (0:8)],
        ytickformat = "{:.1f}",
    )
    for (i, history_i) in history 
        scatter!(ax, 1:length(history_i), history_i, label=L"$\Lambda_%$(i)$ Est.")
        hlines!(ax, [Λ_true[i]], label=L"$\Lambda_%$(i)$ True")
        axislegend(ax)
    end

    save(outfile, fig)
end

"""
Perform the orthogonal iteration on the matrix `A`, for `k_max` iterations.
Prints the diagonal of the matrix `R` at each iteration.
If `Λ_true` and `plot_conv_indices` are nonempty:
    Plots the convergence of the i'th estimated eigenvalue to the i'th true eigenvalue, `Λ_true[i]`.
"""
function orthogonal_iteration(A, k_max, Λ_true = Float64[], plot_conv_indices = Float64[])
    history = Dict(i => Float64[] for i in plot_conv_indices)
    Q = Matrix{Float64}(I, size(A)) # Q starts as the identity matrix.
    for k = 1:k_max
        # Orthogonal iteration:
        Q, R = qr(A * Q)
        # Monitoring:
        Λ_est = diag(R)
        println(round.(Λ_est, sigdigits=4))
        for i in plot_conv_indices
            push!(history[i], Λ_est[i])
        end
    end

    if !isempty(plot_conv_indices) && !isempty(Λ_true)
        plot_convergence(history, Λ_true)
    end 
end

m = 8
Λ_true = 3. .^(0:m-1)
A = random_matrix(m, Λ_true)
orthogonal_iteration(A, 5, reverse(Λ_true), [1, 2, 3])
