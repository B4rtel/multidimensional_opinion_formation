alpha = [[1/2,1/2],[4.0/5,1.0/5]]
alphas="$(join(round.(alpha[1]; digits=2), "-"))_$(join(round.(alpha[2]; digits=2), "-"))"
f0_choice = "stddev_varalp"
initial_density = (x,alp) -> ifelse(alp[1]< 0.6,ifelse(x==[-1,-1], 1.0/30, ifelse(x==[1,-1],14.0/15,0)),ifelse(x==[-5.0/6,1],1.0/30,0))
pa="var2"
function p_a(x, y, alpha, beta, a)
    d = length(x)
    return beta * abs(x[a] - y[a]) + (1 - beta) * sum(alpha[b] * abs(x[b] - y[b]) for b in 1:d)
end
include("teilchenmod_differentdistances.jl")
include("teilchenmod_plots_d_var.jl")


alpha = [[1/2,1/2]]
alphas="$(join(round.(alpha[1]; digits=2), "-"))"
f0_choice = "stddev_fixalp"
initial_density = (x,alp) -> ifelse(x==[-1,-1] || x==[-5.0/6,1], 1.0/30, ifelse(x==[1,-1],14.0/15, 0))
pa="var2"
function p_a(x, y, alpha, beta, a)
    d = length(x)
    return beta * abs(x[a] - y[a]) + (1 - beta) * sum(alpha[b] * abs(x[b] - y[b]) for b in 1:d)
end
include("teilchenmod_differentdistances.jl")
include("teilchenmod_plots_d_var.jl")


alpha = [[4.0/5,1.0/5]]
alphas="$(join(round.(alpha[1]; digits=2), "-"))"
f0_choice = "stddev_fixalp"
initial_density = (x,alp) -> ifelse(x==[-1,-1] || x==[-5.0/6,1], 1.0/30, ifelse(x==[1,-1],14.0/15, 0))
pa="var2"
function p_a(x, y, alpha, beta, a)
    d = length(x)
    return beta * abs(x[a] - y[a]) + (1 - beta) * sum(alpha[b] * abs(x[b] - y[b]) for b in 1:d)
end
include("teilchenmod_differentdistances.jl")
include("teilchenmod_plots_d_var.jl")


alphas="no_alpha"
f0_choice = "stddev_noalp"
initial_density = (x,alp) -> ifelse(x==[-1,-1] || x==[-5.0/6,1], 1.0/30, ifelse(x==[1,-1],14.0/15, 0))
pa="var2"
pa="hk"
function p_a(x, y, alpha, beta, a)
    return  abs(x[a] - y[a])
end
include("teilchenmod_differentdistances.jl")
include("teilchenmod_plots_d_var.jl")


alphas="no_alpha"
f0_choice = "stddev_noalp"
initial_density = (x,alp) -> ifelse(x==[-1,-1] || x==[-5.0/6,1], 1.0/30, ifelse(x==[1,-1],14.0/15, 0))
pa="euklidean"
function p_a(x, y, alpha, beta, a)
    d = length(x)
    return  sqrt(sum(abs(x[b] - y[b])^2 for b in 1:d))
end
include("teilchenmod_differentdistances.jl")
include("teilchenmod_plots_d_var.jl")
