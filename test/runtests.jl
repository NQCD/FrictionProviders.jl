using SafeTestsets

@safetestset "Cube" begin include("test_cube.jl") end
@safetestset "CubeLDFA" begin include("test_cube_density.jl") end
@safetestset "SciKitLDFA" begin include("test_scikit_density.jl") end
