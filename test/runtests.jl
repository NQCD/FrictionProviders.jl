using SafeTestsets

@safetestset "Cube" begin include("test_cube.jl") end
@safetestset "CubeDensity" begin include("test_cube_density.jl") end
@safetestset "SciKitDensity" begin include("test_scikit_density.jl") end
# @safetestset "ACEdsFrictionTensor" begin include("test_aceds_odf.jl") end