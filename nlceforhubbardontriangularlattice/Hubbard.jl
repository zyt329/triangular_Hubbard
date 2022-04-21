using LinearAlgebra
using SparseArrays
using Arpack

"""
    Check for the state, how many spins points up.
    Input: state: the state to check.
           N: the number of sites in the cluster.

return the number of spin-up and spin-down particles in the state: (N_up, N_dn)
"""
function chk_up_dn(state::Int64,N::Int64)
    state_binary = digits!(zeros(Int64, 64), state, base = 2)[1:(2*N)]
    up_state_binary = state_binary[1:2:(2*N)]
    dn_state_binary = state_binary[2:2:(2*N)]
    return (sum(up_state_binary), sum(dn_state_binary))
end

"""
    Loop over all states, check to which sector the state belongs.
return a dictionary:
    :states::Array{Array{Int,1},2} : states[N_up][N_dn] gives an array of states in the (N_up, N_dn) sector
    :state_tot::Array{Int,2} : state_tot[N_up][N_dn] gives the total number of states in the (N_up, N_dn) sector.
    :state_num::Dict{Int64, Int64} : give the state as key, returns the numbering of the state in its (N_up, N_dn) sector.
"""
function sectors_info_gen(;N::Int64)
    states::Array{Array{Int,1},2} = Array{Int,1}[[] for N_up=0:N, N_dn=0:N]
    state_tot::Array{Int,2} = Int[0 for N_up=0:N, N_dn=0:N]
    state_num::Dict{Int64, Int64} = Dict{Int64, Int64}()
    for state in 0:(4^N-1)
        (N_up,N_dn) = chk_up_dn(state, N)
        push!(states[N_up+1, N_dn+1], state)
        state_tot[N_up+1, N_dn+1] += 1
        state_num[state] = state_tot[N_up+1, N_dn+1]
    end
    @assert(sum(state_tot) == 4^N, "total number of states is not 4^N")
    return Dict{Symbol, Any}(:states => states, :state_tot => state_tot, :state_num => state_num)
end

function update_val(row_inds, col_inds, vals;row_ind, col_ind, val)
    push!(row_inds, row_ind)
    push!(col_inds, col_ind)
    push!(vals, val)
end

function H_sector(;t::Real, U::Real, N::Int64, N_up::Int64, N_dn::Int64, sectors_info::Dict{Symbol,Any}, bonds)
    row_inds = Int64[]
    col_inds = Int64[]
    vals = Float64[]
    states = sectors_info[:states][N_up+1, N_dn+1]
    state_tot = sectors_info[:state_tot][N_up+1, N_dn+1]
    state_num = sectors_info[:state_num]
    for state in states #loop over all states in the sector
        state_binary = digits!(zeros(Int64, 64), state, base = 2)[1:(2*N)]
        state_binary_up = state_binary[1:2:(2*N)]
        state_binary_dn = state_binary[2:2:(2*N)]
        #==================== Hubbard U term =====================#
        Hu = 0
        for site in 1:N
            Hu += U*state_binary_up[site]*state_binary_dn[site]
        end
        update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[state], val = Hu)
        # Delete chemical potential in hamiltonian, restore it when calculating
        # thermal average.
        #================= chemical potential term =================
        update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[state], val = -μ * (N_up+N_dn))=#
        #================= Kinetic term ===============#
        for bond in bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]; s2 = bond[2]
            #================== spin-up entries ================#
            if (state_binary_up[s1], state_binary_up[s2])==(1,0) || (state_binary_up[s1], state_binary_up[s2])==(0,1)
                sgn = (-1)^sum(state_binary[(2*s1):(2*(s2-1))])
                flipped_state = state ⊻ (1<<((2*s1-1)-1))
                flipped_state = flipped_state ⊻ (1<<((2*s2-1)-1))
                update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[flipped_state], val = sgn * (-t))
            end
            #================== spin-dn entries ================#
            if (state_binary_dn[s1], state_binary_dn[s2])==(1,0) || (state_binary_dn[s1], state_binary_dn[s2])==(0,1)
                sgn = (-1)^sum(state_binary[(2*s1+1):(2*s2-1)])
                flipped_state = state ⊻ (1<<(2*s1-1))
                flipped_state = flipped_state ⊻ (1<<(2*s2-1))
                update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[flipped_state], val = sgn * (-t))
            end
        end
    end
    return sparse(row_inds, col_inds, vals, state_tot, state_tot, +)
end
