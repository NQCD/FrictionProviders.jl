using FrictionProviders
using Test
using NQCBase
using PythonCall

"""
import sys
sys.path.insert(0, ".")
"""
py_sys = pyimport("sys")
py_sys.path.insert(0, ".")

cube = pyimport("cube")

filename = "test.cube"

c = FrictionProviders.Cube(filename)
cube_object = cube.cube()
cube_object.read(filename)

@testset "Compare to cube.py" begin
    python_density = permutedims(reshape(pyconvert(Vector{Float64}, cube_object."density"), pyconvert(Int, cube_object.x_len), pyconvert(Int, cube_object.y_len), pyconvert(Int, cube_object.z_len)), (3,2,1)) # Permutation as Julia reshape works row-major
    @info size(python_density)
    @info size(c.density)
    @test c.density == python_density

    r = c.origin .+ 1e-8
    @test c(r) ≈ pyconvert(Float64, cube_object(au_to_ang.(r)...))

    for _ = 1:10
        r = c.origin
        r += rand() * c.cell.vectors[:, 1]
        r += rand() * c.cell.vectors[:, 2]
        r += rand() * c.cell.vectors[:, 3]

        @test c(r) ≈ pyconvert(Float64, cube_object(au_to_ang.(r)...))
    end

    r = c.origin + c.cell.vectors[:, 1] + c.cell.vectors[:, 2] + c.cell.vectors[:, 3] .- 1e-8
    @test c(r) ≈ pyconvert(Float64, cube_object(au_to_ang.(r)...))
end
