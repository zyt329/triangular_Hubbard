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
function E_Quants(;T::Float64, N::Int64, sectors_info::Dict{Symbol,Any}, bonds::Array{Array{Int,1},1})
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
        for (ind, state)) in enumerate(states)
            push!(Esq, Es[ind]^2)
            push!(M, ((N_up-N_dn)*1/2))
            push!(Msq, ((N_up-N_dn)*1/2)^2)
            push!(N_tot, N_up+N_dn)
        end

    end
    return (P,E,Esq,M,Msq)
end
