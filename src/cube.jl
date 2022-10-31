
using StaticArrays: SVector

struct Cube{T}
    origin::SVector{3,T}
    shape::Tuple{Int,Int,Int}
    cell::PeriodicCell{T}
    density::Array{T,3}
end

function Cube(filename, ::Type{T}=Float64) where {T<:AbstractFloat}

    data = readdlm(filename; skipstart=2)
    natoms = convert(Int, data[1,1])
    origin = SVector{3,T}(data[1,2:4])
    shape = Tuple(data[2:4])
    vectors_rows = data[2:4,2:4] .* shape
    vectors = permutedims(vectors_rows, (2,1))
    cell = PeriodicCell{T}(vectors, [true, true, true])

    volumetric_data = permutedims(data[5+natoms:end,:], [2,1])
    indices = volumetric_data .!= ""
    density = convert(Vector{T}, volumetric_data[indices])
    density = permutedims(reshape(density, reverse(shape)), 3:-1:1)

    Cube(origin, shape, cell, density)
end

function (cube::Cube)(r::AbstractVector)
    r = r - cube.origin
    r = cube.cell.inverse * r

    indices = floor.(Int, r .* cube.shape) .+ 1

    if any(indices .< 1) || any(indices .> cube.shape)
        throw(DomainError(r, "cannot evaluate the density outside of the cell."))
    end

    return cube.density[indices...]
end
