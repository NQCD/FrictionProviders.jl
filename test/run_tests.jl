using NQCModels: Free, FrictionModels
using PyCall
using NQCBase
using CubeLDFAModel
#using SciKitLDFAModel
using Test
using Unitful
using Pandas


include("test_cube.jl")

aseio = pyimport("ase.io")
ase_atoms = aseio.read("start.in")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

model = LDFAModel(Free(), "test.cube", atoms, cell; friction_atoms=[1, 2])

@testset "friction!" begin
    F = zeros(6, 6)
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(model, F, R)

    for _=1:10
        r .= 0
        r += rand() * cell.vectors[:,1]
        r += rand() * cell.vectors[:,2]
        r += rand() * cell.vectors[:,3]

        FrictionModels.friction!(model, F, R)
    end

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3]
    FrictionModels.friction!(model, F, R)

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3] + rand(3)
    FrictionModels.friction!(model, F, R)
end

# SCI-KIT MODEL
dscr_d = pyimport("dscribe.descriptors")
model_ml = read_pickle("scikit_model_h2cu.sav")
ase_atoms = aseio.read("h2cu_start.in")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

desc = dscr_d.SOAP(
    species = ["Cu", "H"],
    periodic = true,
    rcut = 2.5,
    nmax = 1,
    lmax = 1,
    average="off" 
)
model = LDFAModel(Free(), atoms, cell, desc, model_ml; friction_atoms=[55, 56])

@testset "ScikitModel!" begin
    F = zeros(56, 56)
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(model, F, R)

    for _=1:10
        r .= 0
        r += rand() * cell.vectors[:,1]
        r += rand() * cell.vectors[:,2]
        r += rand() * cell.vectors[:,3]

        FrictionModels.friction!(model, F, R)
    end

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3]
    FrictionModels.friction!(model, F, R)

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3] + rand(3)
    FrictionModels.friction!(model, F, R)
end