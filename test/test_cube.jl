using FrictionProviders
using Test
using NQCBase
using PythonCall


py"""
import sys
sys.path.insert(0, ".")
"""

cube = pyimport("cube")

filename = "test.cube"

c = FrictionProviders.Cube(filename)
cube_object = cube.cube()
cube_object.read(filename)

@testset "Compare to cube.py" begin
    @test c.density == cube_object."density".reshape(cube_object.x_len, cube_object.y_len, cube_object.z_len)

    r = c.origin .+ 1e-8
    @test c(r) ≈ cube_object(au_to_ang.(r)...)

    for _=1:10
        r = c.origin
        r += rand() * c.cell.vectors[:,1]
        r += rand() * c.cell.vectors[:,2]
        r += rand() * c.cell.vectors[:,3]

        @test c(r) ≈ cube_object(au_to_ang.(r)...)
    end

    r = c.origin + c.cell.vectors[:,1] + c.cell.vectors[:,2] + c.cell.vectors[:,3] .- 1e-8
    @test c(r) ≈ cube_object(au_to_ang.(r)...)
end

