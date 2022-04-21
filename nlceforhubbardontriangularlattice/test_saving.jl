using JLD
using DelimitedFiles

quant_names = ['a']
sample = [[]]
for i in 1:2^18
    push!(sample[1], i)
end

"""
Quantities[n][i] is the nth quantities measured at Temperature i.
"""
function printing(Quantities; Quant_names = quant_names, name = "test")
    println(Quantities)
    open(
        "E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/"*name*".dat",
        "w",
    ) do io
        writedlm(io, quant_names)
        for (i,Quant) in enumerate(Quantities[1])
            writedlm(io, [Quantities[n][i] for n in 1:length(Quantities)]')
        end
    end

end

printing(sample; name = "test_print2txt")
