include("Functions.jl");
using BenchmarkTools
using Plots
using StatsPlots
M = Model(0.1, 0.1, 100, 100, -100., 1.);
G = CreateGrid(M, "Cold")


Initial_Energy = V(M, G) + T(M, G)
Initial_Action = M.dt * sum(Initial_Energy)

Total_Action = Initial_Action;

@inbounds for i in 1:3000
    @inbounds for j in 1:M.Nt-1
        @inbounds for k in 1:M.Nx-1

            G, Total_Action = Iterate(j, k, 0.5, G, M, Total_Action)
            if Total_Action != Total_Action
                println("NaN detected.")
            end
        end
    end
end

TotalEnergy = zeros(Float64, M.Nt, M.Nx)
TotalField = zeros(Float64, M.Nt, M.Nx)

A = @animate for i in 1:50000
    @inbounds for j in 1:M.Nt-1
        @inbounds for k in 1:M.Nx-1

            G, Total_Action = Iterate(j, k, 0.5, G, M, Total_Action)
            if Total_Action != Total_Action
                println("NaN detected.")
            end
        end
    end
    Energy = V(M, G)
    TotalField += G
    TotalEnergy += Energy
    v = vec(Energy)
    phi = vec(G)
    #println(minimum(v)/i, " ", maximum(v)/i)
    w = histogram(phi/i)
    #w = Plots.plot!([sqrt(-M.m/M.λ)], seriestype = :vline, label = "Vacuum Expectation Energy")
    #println(mean(phi/i))
    l = heatmap(G, ylabel = "Grid Time", xlabel = "Grid Position", colorbar = true, legend = false, title = "Field Values of Scalar 1+1D Field Theory \n for Markov Run $(i), m^2 = $(M.m), λ = $(M.λ)")
    g = density(phi, bins = 200, xlabel = "Field Value", title = "Distribution of Field Values", colorbar = true, legend = false)
    plot(l, g, layout = (2, 1), size = (700, 700))
    #plot(l, title = "Field Values of Scalar 1+1D Field Theory \n for Markov Run $(i), m^2 = $(M.m), λ = $(M.λ)")
end every 20
gif(A, "anim_fps15.gif", fps = 100)

heatmap(G .* G)
Initial_Energy = V(M, G) + T(M, G)
heatmap(Initial_Energy)
