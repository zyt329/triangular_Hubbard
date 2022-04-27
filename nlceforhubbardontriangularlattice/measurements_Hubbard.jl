#= functions to
    1. Diagonalize hamiltonian for the sectors of each cluster
    2. For each eigenstate, calculate expectation values of desired quantities
    3. Store the eigenvalues, <quantities> of each cluster in a file
    4. For each parameter(T,μ,h...), calculate Thermal average of desired quantities and feed back to the NLCE Python script.
=#
"""
    Function to produce eigenvalues(energies) of a cluster and corresponding expectation values of quantities of interest.
    Input:
        N::Float64 : number of sites in the cluster.
        sectors_info::Dict{Symbol,Any} : information of symmetry sectors of the cluster.
        bonds::Array{Array{Int,1},1} : bond information of the cluster.
    Output:
        [E, Quantity1, Quantity2...]. E, Quantity1, Quantity2... are all arrays of length 2^N
"""
function E_Quants(;N::Int64, U::Real, t::Real, sectors_info::Dict{Symbol,Any}, bonds::Array{Array{Int,1},1})
    #P::Float64 = 0
    E = Float64[]
    Esq = Float64[]
    M = Float64[]
    Msq = Float64[]
    N_tot = Float64[]
    for N_up=0:N, N_dn=0:N
        #println("N_up=$N_up, N_dn=$N_dn")
        # construct hamiltonian for the sector and diagonalize
        H_m = H_sector(;t=t, U=U, N=N, N_up=N_up, N_dn=N_dn, sectors_info=sectors_info, bonds=bonds)
        (Es, states) = eigen(Hermitian(Matrix(H_m)))
        # calculate quantiies for each eigenstate
        for (ind, En) in enumerate(Es)
            push!(E, En)
            push!(Esq, En^2)
            push!(M, (N_up-N_dn)*1/2)
            push!(Msq, ((N_up-N_dn)*1/2)^2)
            push!(N_tot, (N_up+N_dn))
        end
    end
    return [E, Esq, M, Msq, N_tot]
end



using JLD
using DelimitedFiles


#=sample = [Float64[], Float64[], Float64[], Float64[], Float64[]]
for i in 1:5
    push!(sample[i], i)
    push!(sample[i], 2i)
end=#


"""
    prints the Quantities to a .dat file for save in a given folder for each cluster. Saved name is name_Ntop
    Input:
        Quantities::Array{Array{Float64,1},1} : An array of Arrays of quantities
            example: sample = [Float64[], Float64[], Float64[],     Float64[], Float64[]]
                for i in 1:5
                    push!(sample[i], i)
                    push!(sample[i], 2i)
                end
        Quant_names::Array{Any,2} : An array of names of quantities
            example: ['E' "Esq" 'M' "Msq" "N_tot"]
        name::string : name of the saved files
        NTOP::Int64 : topological number of the cluster whose quantities we're saving.
"""
function printing(Quantities; Quant_names = quant_names, name = "test", NTOP::Int64,N::Int64)
    #println(Quantities)
    open(
        "./Hubbard_data/"*name*"_$(N)_$(NTOP)"*".dat",
        "w",
    ) do io
        writedlm(io, quant_names)
        for (i,Quant) in enumerate(Quantities[1])
            writedlm(io, [Quantities[n][i] for n in 1:length(Quantities)]')
        end
    end
end


# Quantities = E_Quants(;N=N, sectors_info=sectors_info, bonds=bonds)

# printing(Quantities; name = "test_print2txt", NTOP=1)


"""
    do thermal average of quantities read from the file created by E_Quants() and printing()
    return:
        An array of thermal average of quantities at temperature of T
"""
function thermal_avg(;T::Real, μ::Real, name::String = "test", NTOP::Int64, N::Int64)
    # read in the data from file
    global quantities
    open(
        "./Hubbard_data/"*name*"_$(N)_$(NTOP)"*".dat",
        "r",
    ) do io
        quantities=readdlm(io, skipstart=1) # skip the first line
    end
    E::Array{Float64, 1} = quantities[:, 1]
    Esq::Array{Float64, 1} = quantities[:, 2]
    M::Array{Float64, 1} = quantities[:, 3]
    Msq::Array{Float64, 1} = quantities[:, 4]
    N_tot::Array{Float64, 1} = quantities[:, 5]
    # calculate thermal average
    β = 1/T
    Z::Float64 = 0
    E_avg::Float64 = 0
    Esq_avg::Float64 = 0
    M_avg::Float64 = 0
    Msq_avg::Float64 = 0
    N_tot_avg::Float64 = 0
    for (n, En) in enumerate(E)
        P = exp(-β*(En-μ*N_tot[n]))
        Z += P
        E_avg += En * P
        Esq_avg += Esq[n] * P
        M_avg += M[n] * P
        Msq_avg += Msq[n] * P
        N_tot_avg += N_tot[n] * P
    end
    return [Z E_avg Esq_avg M_avg Msq_avg N_tot_avg]
end

quant_names = ['E' "Esq" 'M' "Msq" "N_tot"]
#=
Main.include("Hubbard.jl")
bonds=[[1,2]]
N=2; NTOP=0
sectors_info = sectors_info_gen(N=N)
quantities = E_Quants(N=N, U=10.0, t=1.0, sectors_info=sectors_info, bonds=bonds)
#printing(quantities; Quant_names = quant_names, name = "test", NTOP=NTOP)
thermal_avg(;T=100,μ=5.0, name = "test",N=N, NTOP=NTOP)=#
