using Plots
using Printf

default(titlefont=font(15), guidefont=font(14), tickfont=font(12), legendfont=font(12))

Plots.gr_cbar_width[]=0.0075

grad = cgrad(:sunset, scale=(x -> x^3.5), rev=true)  # x^0.1 exaggerates the low end
grad_s = "x^3.5"

xlims_sol = (-1, 1)
ylims_sol = (-1, 1)

# bin particles
function bin_particles(positions, weights; x_grid, tol=1e-12)
    N = length(x_grid)
    d = size(positions, 1)
    # Create a d-dimensional array for binned weights
    binned_weights = zeros(ntuple(_ -> N, d))
    # Bin each particle
    for i in 1:size(positions, 2)
        inds = ntuple(a -> clamp(findmin(abs.(x_grid .- positions[a, i]))[2], 1, N), d)
        binned_weights[inds...] += weights[i]
    end
    # Collect bin centers and weights for nonzero bins
    idxs = CartesianIndices(binned_weights)
    centers = [Float64[] for _ in 1:d]
    ws = Float64[]
    for I in idxs
        w = binned_weights[I]
        if w > tol
            for a in 1:d
                push!(centers[a], x_grid[I[a]])
            end
            push!(ws, w)
        end
    end
    # Return a tuple of vectors for each dimension and the weights
    return (centers..., ws)
end

positions_init = sol.u[1]
positions_final = sol.u[end]
weights_init = fill(1.0 / size(positions_init, 2), size(positions_init, 2))
weights_final = fill(1.0 / size(positions_final, 2), size(positions_final, 2))

result_init = bin_particles(positions_init, weights_init; x_grid=x_grid)
result_final = bin_particles(positions_final, weights_final; x_grid=x_grid)

centers_init = result_init[1:end-1]
ws_init = result_init[end]
centers_final = result_final[1:end-1]
ws_final = result_final[end]

clim = maximum([maximum(ws_init), maximum(ws_final)])
clim_min = 0

p3 = heatmap(1:2,1:2,[0 0;0 0], colorbar=true,clim = (clim_min,clim), color = grad, grid=false, axis=false, size=(10,500),right_margin=-10Plots.mm,left_margin=-10Plots.mm)
# arrange marginals horizontally with programmatic widths
frc = 7.9  # total figure width factor (d main panels + 1 colorbar panel)
main_width_1 = (frc-1) / (frc * 2)   # relative width for each main panel
cb_width = 1 / frc                 # relative width for the colorbar panel
widths = vcat(fill(main_width_1, 2), cb_width)   # length == d+1


if d == 2
    xs_init, ys_init = centers_init
    xs_final, ys_final = centers_final
    p1 = scatter(xs_init, ys_init, marker_z=ws_init, color=grad, clim=(clim_min, clim),
        xlabel="x1", ylabel="x2", title="Initial particle density", markersize=13, marker=:rect,
        aspect_ratio=:equal, size=(600, 600), colorbar=false, legend=false,
        markerstrokewidth=0, xlims=xlims_sol, ylims=ylims_sol,
        left_margin=7Plots.mm,  tick_direction=:out,widen=true)
    p2 = scatter(xs_final, ys_final, marker_z=ws_final, color=grad, clim=(clim_min, clim),
        xlabel="x1", ylabel="x2", title="Final particle density", markersize=13, marker=:rect,
        aspect_ratio=:equal, size=(600, 600), colorbar=false, legend=false,
        markerstrokewidth=0, xlims=xlims_sol, ylims=ylims_sol,left_margin=2Plots.mm,tick_direction=:out,  widen=true)

    plot(p1, p2, p3, layout = grid(1, 3,  widths=tuple(widths...)), size=(1200+80, 600))
    savefig("plots/initial_final_f_pa$(pa)_d$(d)_N$(N)_beta$(beta)_alphas_$(alphas)_phi$(ϕ_choice)_T$(T)_f0$(f0_choice)_algo$(algo_s)_density_scale_$(density_scale)_N_particles$(N_particles)_scaling$(grad_s).pdf")

elseif d == 3
    xs_init, ys_init, zs_init = centers_init
    xs_final, ys_final, zs_final = centers_final
    p1 = scatter3d(xs_init, ys_init, zs_init, marker_z=ws_init, color=grad, clim=(clim_min, clim),
        xlims=xlims_sol, ylims=ylims_sol, zlims=xlims_sol,
        xlabel="x1", ylabel="x2", zlabel="x3", title="Initial particle density", markersize=9,
        size=(600, 600), colorbar=false, legend=false, markerstrokewidth=0, left_margin=7Plots.mm, tick_direction=:out)
    p2 = scatter3d(xs_final, ys_final, zs_final, marker_z=ws_final, color=grad, clim=(clim_min, clim),
        xlims=xlims_sol, ylims=ylims_sol, zlims=xlims_sol,
        xlabel="x1", ylabel="x2", zlabel="x3", title="Final particle density", markersize=9,
        size=(600, 600), colorbar=false, legend=false, markerstrokewidth=0, tick_direction=:out)
    #plot(p1, p2, layout=(1,2), size=(1200, 600))
    plot(p1, p2, p3, layout = grid(1, 3,  widths=tuple(widths...)), size=(1200+80, 600))
    savefig("plots/initial_final_f_pa$(pa)_d$(d)_N$(N)_beta$(beta)_alphas$(alphas)_phi$(ϕ_choice)_T$(T)_f0$(f0_choice)_algo$(algo_s)_density_scale_$(density_scale)_N_particles$(N_particles)_scaling$(grad_s).pdf")

end

# marginal plots
num_timesteps = length(sol.t)
full_time = sol.t
f_marginals = [zeros(length(x_grid), num_timesteps) for a in 1:d]

for t in 1:num_timesteps
    positions = sol.u[t]
    weights = fill(1.0 / size(positions, 2), size(positions, 2))
    result = bin_particles(positions, weights; x_grid=x_grid)
    centers = result[1:end-1]
    ws = result[end]
    for a in 1:d
        # Bin index for each dimension
        for k in 1:length(ws)
            ia = clamp(findmin(abs.(x_grid .- centers[a][k]))[2], 1, length(x_grid))
            f_marginals[a][ia, t] += ws[k]
        end
    end
end

minval = 0
maxval = maximum([maximum(f_marginals[a]) for a in 1:d])

plot_list = []
for a in 1:d
    push!(plot_list, heatmap(
        x_grid,
        full_time,
        f_marginals[a]' .+ minval,
        xlabel = "x$a", ylabel = "time", title = "Marginal in dimension $a",
        interpolation = :none,
        size = (900, 500),
        #colorbar = true,
        colorbar = false,
        color = grad,
        left_margin=7Plots.mm,
        clim = (minval, maxval),
        right_margin=0Plots.mm,
        bottom_margin=7Plots.mm,
        tick_direction=:out,
        ylims=(0, T)
    ))
end
# add colorbar for a=d
    push!(plot_list, 
     heatmap(1:2,1:2,[0 0;0 0], colorbar=true,clim = (minval, maxval), color = grad, grid=false, axis=false, size=(10,500),right_margin=-10Plots.mm,left_margin=-10Plots.mm))



# arrange marginals horizontally with programmatic widths
frc = 7.4  # total figure width factor (d main panels + 1 colorbar panel)
main_width = (frc-1) / (frc * d)   # relative width for each main panel
cb_width = 1 / frc                 # relative width for the colorbar panel
widths = vcat(fill(main_width, d), cb_width)   # length == d+1

# arrange marginals horizontally (1 row, d columns)
dheatmaps = plot(plot_list..., layout = grid(1, d+1,  widths=tuple(widths...)), size = (500*d+100, 500))#, right_margin=15Plots.mm)
savefig(dheatmaps, "plots/dheatmaps_horiz_pa$(pa)_d$(d)_N$(N)_beta$(beta)_alphas$(alphas)_phi$(ϕ_choice)_T$(T)_f0$(f0_choice)_algo$(algo_s)_density_scale_$(density_scale)_N_particles$(N_particles)_horiz.pdf")
 
