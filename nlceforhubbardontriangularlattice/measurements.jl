function measure(;T::Float64, N::Int64, sectors_info::Dict{Symbol,Any}, bonds::Array{Array{Int,1},1})
    β = 1/T
    P::Float64 = 0
    E::Float64 = 0
    Esq::Float64 = 0
    M::Float64 = 0
    Msq::Float64 = 0
    for m in 0:N
        H_m = H_sector(J=1.0, N=N, m=m, sectors_info=sectors_info, bonds=bonds)
        (Es, states) = eigen(Matrix(H_m))
        for (n, En) in enumerate(Es)
            ΔP = exp(-β*En)
            P += ΔP
            E += En * ΔP
            Esq += En^2 * ΔP
            M += ((2m-N)*1/2) * ΔP
            Msq += ((2m-N)*1/2)^2 * ΔP
        end
    end
    return (P,E,Esq,M,Msq)
end
