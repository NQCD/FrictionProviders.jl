using SafeTestsets

@safetestset "Cube" begin include("test_cube.jl") end
@safetestset "CubeLDFA" begin include("test_cube_ldfa.jl") end
@safetestset "SciKitLDFA" begin include("test_scikit_ldfa.jl") end
