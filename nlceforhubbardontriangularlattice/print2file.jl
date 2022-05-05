using DelimitedFiles


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
function print2file(Quantities; Quant_names = quant_names, name = "test",order)
    #println(Quantities)
    open(
        "./final_data/"*name*"_U=$(U)_order=$(order)"*".dat",
        "w",
    ) do io
        writedlm(io, quant_names)
        writedlm(io, Quantities)
        #=
        for (i,Quant) in enumerate(Quantities[1])
            #writedlm(io, [Quantities[n][i] for n in 1:length(Quantities)]')
        end=#
    end
end

quant_names = ["Temp" 'E' 'M' "C" "Chi" "N_tot" "lnZ" "S"]

for N in order2print
    # preprocess raw_O to make the 3rd and 4th quantities C and Chi
    raw_O[3,N+1,:] = raw_O[3,N+1,:] ./ Temps.^2
    raw_O[4,N+1,:] = raw_O[4,N+1,:] ./ Temps
    quant2print = hcat(Temps, raw_O[:,N+1,:]')
    # add entropy S to the last quantity
    Ss = raw_O[6,N+1,:] + 1 ./ Temps .* (raw_O[1,N+1,:]- μ * raw_O[5,N+1,:])
    quant2print = hcat(quant2print, Ss)
    print2file(quant2print; Quant_names = quant_names, name = "tri_Hub", order = N)
end
