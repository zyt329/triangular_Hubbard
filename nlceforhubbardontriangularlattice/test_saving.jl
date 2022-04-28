using JLD
using DelimitedFiles

quant_names = ['a' 'b' 'c' 'd' 'e']
sample = [[],[],[],[],[]]
for i in 1:2^2
    push!(sample[1], rand())
    push!(sample[2], rand())
    push!(sample[3], rand())
    push!(sample[4], rand())
    push!(sample[5], rand())
end

"""
Quantities[n][i] is the nth quantities measured at Temperature i.
"""
function printing(Quantities; Quant_names = quant_names, name = "test")
    println(Quantities)
    open(
        "./"*name*".dat",
        "w",
    ) do io
        writedlm(io, quant_names)
        for (i,Quant) in enumerate(Quantities[1])
            writedlm(io, [Quantities[n][i] for n in 1:length(Quantities)]')
        end
    end

end

printing(sample; name = "test_print2txt1")

@time open(
    "./test_print2txt1"*".dat",
    "r",
) do io
    quantities=readdlm(io, skipstart=1) # skip the first line
end

function reading_quantities(;name::String = "order=7_U=12", NTOP::Int64, N::Int64)
    global quantities
    @time open(
        "E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/Hubbard_data/"*name*"_$(N)_$(NTOP)"*".dat",
        "r",
    ) do io
        quantities=readdlm(io, skipstart=1) # skip the first line
    end
    return quantities
end

quantities = reading_quantities(name= "order=7_U=12", NTOP=0, N=7)

"""
    do thermal average of quantities read from the file created by E_Quants() and printing()
    return:
        An array of thermal average of quantities at temperature of T
"""
function thermal_avg(;T::Real, μ::Real, quantities)
    # passing in quantities
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

function thermal_avg_loop(;Temps, μ::Real, quantities)
    Zs::Array{Float64,1} = 0
    E_avgs::Array{Float64,1} = 0
    Esq_avgs::Array{Float64,1} = 0
    M_avgs::Array{Float64,1} = 0
    Msq_avgs::Array{Float64,1} = 0
    N_tot_avgs::Array{Float64,1} = 0
    for T in Temps
        thermal_avg(;T=T, μ=μ, quantities = quantities)
    end
    return [Z, E_avgs, Esq_avgs, M_avgs, Msq_avgs, N_tot_avgs]
end




@time for T in range(0.1, 10, length = 100)
    thermal_avg(;T=T,μ=5.0, quantities = quantities)
end
