import Pkg.instantiate
instantiate()

using LinearAlgebra
using CairoMakie

"""
Create a diagonalizable matrix of size m x m, with the provided eigenvalues `Λ`. 
Performs a QR decomposition of a random matrix to create an orthogonal matrix `Q`,
and returns `Q diag(Λ) Q`.
Arguments:
    * `m` is the size of the matrix.
    * `Λ` is a vector of the eigenvalues.
"""
function random_matrix(m, Λ)
    Q, _ = qr(randn(m, m))
    return Q * Diagonal(Λ) * Q'
end

"""
Plot the convergence of the estimated eigenvalues to the true eigenvalues.
The plot has two subplots:
  * The first subplot shows the estimated and true eigenvalues.
  * The second subplot shows the actual and theoretical convergence rates.
Arguments:
  * `Λ_history` is a dictionary mapping the index of the eigenvalue to a vector of the estimated eigenvalues at each iteration.
  * `Λ_true` is a vector of the true eigenvalues.
  * `outfile` is the path of the file to save the plot to.
"""
function plot_convergence(Λ_history::Dict{Int64, Vector{Float64}}, Λ_true::Vector{Float64}, outfile="convergence.png")
    fig = Figure()
    indices = Λ_history |> keys |> collect |> sort
    indices_str = join(indices, ", ")
    values_ax = Axis(fig[1, 1],
        title = L"Estimated vs. true $\Lambda_i$ for $p \in [%$indices_str]",
        xlabel = L"Iteration $k$",
        ylabel = "Estimated eigenvalue",
        yscale = Makie.Symlog10(3.0), # Symlog scale allows for negative values.
        yticks = [-3 .^ (8:-1:1); 0; 3 .^ (1:8)],
        ytickformat = "{:d}",
    )
    conv_ax = Axis(fig[2, 1],
        title = L"Actual vs. theoretical convergence of $\Lambda_p$ for $p \in [%$indices_str]",
        xlabel = L"Iteration $k$",
        ylabel = "Convergence",
        yscale = log10,
    )

    linestyles = [:solid, :dash]
    colors = map(i ->  Makie.wong_colors()[i], 1:length(indices))
    m = length(Λ_true)
    for (i, (p, p_history)) in enumerate(Λ_history)
        n = length(p_history)
        color = colors[i]

        # Plot the estimated and true eigenvalues.
        lines!(values_ax, 1:n, p_history, linestyle=linestyles[1], color=color)
        hlines!(values_ax, [Λ_true[p]], linestyle=linestyles[2], color=color)
 
        # Plot the actual and theoretical convergence rates.
        conv = abs.(p_history .- Λ_true[p]) ./ Λ_true[p]
        theoretical_conv_base =
            p == 1 ? abs(Λ_true[2] / Λ_true[1]) :
            p == m ? abs(Λ_true[m] / Λ_true[m - 1]) :
            max(abs(Λ_true[p+1] / Λ_true[p]), abs(Λ_true[p] / Λ_true[p-1]))
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
            ["Estimated", "True"],
            map(p -> L"\lambda_{%$p}", indices),
        ],
        ["Series", "Eigenvalue"]
    )
    Legend(fig[2, 2],
        [group_style, group_color],
        [
            ["Actual", "Theoretical"],
            map(p -> L"\lambda_{%$p}", indices),
        ],
        ["Series", "Eigenvalue"]
    )

    save(outfile, fig)
end

"""
Plot the provided history of 2-norms, along with the analytical estimates for the decay at each iteration.
"""
function plot_norms(norm_history::Dict{Int64, Vector{Float64}}, Λ_true::Vector{Float64}, outfile="norms.png")
    fig = Figure()
    indices = norm_history |> keys |> collect |> sort
    indices_str = join(indices, ", ")
    ax = Axis(fig[1, 1],
        title = L"Actual vs. analytical 2-norms of subblock $(p+1:m,1:p)$ for $p \in [%$indices_str]",
        xlabel = L"Iteration $k$",
        ylabel = L"2-norm of sublock of $A_k = Q_k^T A Q_k$",
        yscale = log10,
    )

    linestyles = [:solid, :dash]
    colors = map(i -> Makie.wong_colors()[i], 1:length(indices))
    m = length(Λ_true)
    for (i, (p, p_history)) in enumerate(norm_history)
        n = length(p_history)
        color = colors[i]

        # Plot the actual and analytical 2-norms.
        norms = p_history
        analytical_est_base =
            p == 1 ? abs(Λ_true[2] / Λ_true[1]) :
            p == m ? abs(Λ_true[m] / Λ_true[m - 1]) :
            abs(Λ_true[p+1] / Λ_true[p])
        analytical_est = analytical_est_base .^ (1:n)
        lines!(ax, 1:n, norms, linestyle=linestyles[1], color=color)
        lines!(ax, 1:n, analytical_est, linestyle=linestyles[2], color=color)
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
            ["Actual", "Analytical"],
            map(p -> L"%$p", indices),
        ],
        ["Series", L"$p$"]
    )

    save(outfile, fig)
end

"""
Perform the orthogonal iteration on the matrix `A`, for `k_max` iterations.
Prints the diagonal of the matrix `R` at each iteration.
If `Λ_true` and `plot_conv_indices` are nonempty:
    Plots the convergence of the eigenvalues for the indices in `plot_conv_indices`.
If `plot_norm_subblock_indices` is nonempty:
    Plots the 2-norm of the subblock (p+1:m,1:p) of the matrices `A_k = Q_k^T A Q_k` for the indices in `plot_norm_subblock_indices`.
"""
function orthogonal_iteration(A, k_max, Λ_true = Float64[], plot_conv_indices = Float64[], plot_norm_subblock_indices = Float64[])
    Λ_history = Dict(p => Float64[] for p in plot_conv_indices)
    norm_history = Dict(p => Float64[] for p in plot_norm_subblock_indices)

    m = size(A, 1)
    Q = Matrix{Float64}(I, (m, m)) # Q starts as the identity matrix.
    for k = 1:k_max
        # Orthogonal iteration:
        Q, R = qr(A * Q)

        # Monitoring:
        Λ_est = diag(R)
        println(round.(Λ_est, sigdigits=4))
        for p in plot_conv_indices
            push!(Λ_history[p], Λ_est[p])
        end

        # Track the 2-norm of each p+1:m,1:p subblock of $A_k = Q_k^T A Q_k$.
        for p in plot_norm_subblock_indices
            A_k = Q' * A * Q
            push!(norm_history[p], norm(A_k[p+1:m, 1:p], 2))
        end
    end

    if !isempty(plot_conv_indices) && !isempty(Λ_true)
        plot_convergence(Λ_history, Λ_true)
    end

    if !isempty(norm_history)
        plot_norms(norm_history, Λ_true)
    end
end

m = 8
Λ_true = 3. .^(0:m-1)
A = random_matrix(m, Λ_true)
orthogonal_iteration(A, 10, reverse(Λ_true), [1, 2, 3], [4])
