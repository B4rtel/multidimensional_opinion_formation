using LinearAlgebra
using Plots
using DifferentialEquations
using ElasticArrays

# parameters
d = 2 # dimension of the opinion space
N = 2^5*3+1
dx = 2.0 / (N-1)
x_grid = collect(-1.0:dx:1.0)
beta = 0.5
T = 2500.0

# ODE solver
algo = Vern9() 
algo_s = "Vern9()"


L = length(alpha)


# Sampling: number of particles at each location proportional to density (rounded)
density_scale = 25
peak_opinions = Vector{Vector{Float64}}()
peak_alphas = Vector{Vector{Float64}}()


for l in 1:L
    # Create a product iterator over x_grid, repeated d times
    for pos_tuple in Iterators.product(ntuple(_ -> x_grid, d)...)
        pos = collect(pos_tuple)  # Convert tuple to vector
        dens = initial_density(pos, alpha[l])
        n_particles = round(Int, dens * density_scale)
        for _ in 1:n_particles
            push!(peak_opinions, pos)
            push!(peak_alphas, alpha[l])
        end
    end
end

N_particles = length(peak_opinions)
println("Total particles: $N_particles")
# All particles have equal weight
peak_weights = fill(1.0 / N_particles, N_particles)

# Interaction function
r_1 = 2.0/5
r_2 = 1.0/2
add = r_2-r_1
q(r) = r^2/(r^2+(1-r)^2)
ϕ(r) = ifelse(r <= r_1, 1.0, ifelse(r < r_2, q((r_2-r)/(r_2-r_1)), 0.0))
ϕ_choice = "smooth_cutoff_r1$(r_1)_add$(add)"


# ODE system setup
function make_ode_sys(peak_opinions, peak_weights, peak_alphas, ϕ, p_a, beta, T)
    p = (;beta, peak_alphas, ϕ, peak_weights)
    u_0 = ElasticArray{Float64}(hcat(peak_opinions...))
    tspan = (0.0, T)
    function odefunc!(du::ElasticArray{Float64}, u::ElasticArray{Float64}, p, t::Float64)
        N_peaks = size(u, 2)
        d = size(u, 1)
        fill!(du, 0.0)
        for i in 1:N_peaks
            for j in 1:N_peaks
                if i != j
                    dist = [p_a(u[:,i], u[:,j], p.peak_alphas[i], beta, a) for a in 1:d]
                    du[:, i] .+= ϕ.(dist) .* (u[:, j] .- u[:, i]) .* p.peak_weights[j]
                end
            end
        end
        du ./= N_peaks
    end
    return odefunc!, u_0, tspan, p
end

# Solve ODE
odefunc!, u_0, tspan, p = make_ode_sys(peak_opinions, peak_weights, peak_alphas, ϕ, p_a, beta, T)
prob = ODEProblem(odefunc!, u_0, tspan, p)
sol = solve(prob, algo, saveat=range(0, T, length=10000))
