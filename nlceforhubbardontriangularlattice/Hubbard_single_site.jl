# calculate properties for single-site Hubbard model
# Should import quantities form python script

# calculate quantities
function single_site_quantities(;Ts, μ::Real, U::Real)
    Z = Float64[]
    E = Float64[]
    Esq = Float64[]
    M = Float64[]
    Msq = Float64[]
    N_tot = Float64[]
    for T in Ts
        β = 1/T
        Zt = 1+2*exp(β*μ)+exp(-β*U+2β*μ)
        push!(Z, Zt)
        push!(E,( U*exp(-β*U+2β*μ) ) / Zt)
        push!(Esq,(U^2*exp(-β*U+2β*μ) ) / Zt)
        push!(M,0)
        push!(Msq , (1/2*exp(β*μ))/Zt)
        push!(N_tot , (2*exp(β*μ)+2*exp(-β*U+2β*μ))/Zt)
    end
    return [Z, E,Esq, M, Msq, N_tot]
end

result = single_site_quantities(Ts=[1.0,100],μ=5.0,U=10.0)
println(result[1])
#=
using Plots

E = single_site_quantities(Ts=range(0.1,10,length=100),μ=0.0,U=10.0)[1]

plot(range(0.1,10,length=100), E)=#
