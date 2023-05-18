using SafeTestsets

@safetestset "Cube" begin include("test_cube.jl") end
@safetestset "CubeDensity" begin include("test_cube_density.jl") end
@safetestset "SciKitDensity" begin include("test_scikit_density.jl") end
@safetestset "SchNetODF" begin include("test_schnet_odf.jl") end