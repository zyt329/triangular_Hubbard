# for square lattice
# c means creation operator, d means annihilation operator

function index(ix::Int64, iy::Int64, L::Int64)
    ind = (ix-1) * L + iy
    return ind
end

function coord(index::Int64, L::Int64)
    coordinate = (Int(floor(index/L))+1 ,mod1(index, L))
    return coordinate
end

"""
    Calculate the index of a site's x-direction neighbor given the site's index
"""
function x_neib(ind::Int64, L::Int64)
    (ix, iy) = coord(ind, L)
    ix_new = mod1(ix+1, L)
    ind_new = index(ix_new, iy, L)
    return ind_new
end

energy(kx, ky, t) = -2t * (cos(kx) + cos(ky))

n_F(E, T, μ) = 1 / (exp((E-μ)/T) + 1)

function k_points(L)
    k_list = []
    for ix in 0:(L-1), iy in 0:(L-1)
        push!(k_list, (-pi+2pi*ix/L, -pi+2pi*iy/L))
    end
    return k_list
end
#println(k_pts(4))

function cidj(;ix::Int64, iy::Int64, jx::Int64, jy::Int64, L::Int64, t::Real, T::Float64, μ::Float64)
    result = 0.0 + im*0.0
    # generate all the k points
    k_pts = k_points(L)
    for (kx, ky) in k_pts
        result += n_F(energy(kx, ky, t), T, μ) * exp(im*( kx*(ix-jx) + ky*(iy-jy) ))
    end
    result = result / L^2
    return result
end
# test_result = cidj(ix=5, iy=5, jx=5, jy=5, L=4, t=1, T=0.01, μ=4.1)


function cidj_list_gen(;L::Int64, t::Real, T::Real, μ::Real)
    cidj_list = Array{Complex{Float64}, 2}(undef, L^2, L^2)
    for i_ind in 1:L^2, j_ind in 1:L^2
        (ix, iy) = coord(i_ind, L)
        (jx, jy) = coord(j_ind, L)
        cidj_list[i_ind, j_ind] = cidj(ix=ix, iy=iy, jx=jx, jy=jy, L=L, t=t, T=T, μ=μ)
    end
    return cidj_list
end

# println(cidj_list_gen(;L=2, t=1, T=0.01, μ=4.1))
function jxijxj(;i_ind::Int64, j_ind::Int64, L::Int64, t::Real, T::Real, μ::Real)
    cidj = cidj_list_gen(L=L, t=t, T=T, μ=μ)
    result = 0.0 + 0.0*im

    result += (i_ind==x_neib(j_ind, L))*cidj[x_neib(i_ind, L), j_ind]
    result -= (x_neib(i_ind, L)==j_ind)*cidj[x_neib(j_ind, L), i_ind]

    result += (i_ind==j_ind)*cidj[x_neib(i_ind, L), x_neib(j_ind, L)]
    result -= (x_neib(i_ind, L)==x_neib(j_ind, L))*cidj[j_ind, i_ind]

    result += (x_neib(i_ind, L)==x_neib(j_ind, L))*cidj[i_ind, j_ind]
    result -= (i_ind==j_ind)*cidj[x_neib(j_ind, L), x_neib(i_ind, L)]

    result += (x_neib(i_ind, L)==j_ind)*cidj[i_ind, x_neib(j_ind, L)]
    result -= (i_ind==x_neib(j_ind, L))*cidj[j_ind, x_neib(i_ind, L)]

    return result
end

function jxijxj_list_gen(;L::Int64, t::Real, T::Real, μ::Real)
    jxijxj_list = Array{Complex{Float64}, 2}(undef, L^2, L^2)
    for i_ind in 1:L^2, j_ind in 1:L^2
        jxijxj_list[i_ind, j_ind] = jxijxj(i_ind=i_ind,j_ind=j_ind, L=L, t=t, T=T, μ=μ)
    end
    return jxijxj_list
end

jxijxj_list = jxijxj_list_gen(;L=4, t=1, T=1.0, μ=2.1)
println(jxijxj_list)
