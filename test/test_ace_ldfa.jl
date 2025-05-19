using Test
using FrictionProviders
using PythonCall
using NQCBase: NQCBase
using NQCModels: FrictionModels
import Random
using Unitful: @u_str
using JuLIP
using ACEpotentials

function ace_model(model_path, cur_atoms)
    IP = read_dict(load_dict(model_path)["IP"])
    JuLIP.set_calculator!(cur_atoms, IP)
    
    model = AdiabaticModels.JuLIPModel(cur_atoms)
    return model
end

# ACE MODEL
aseio = pyimport("ase.io")

ase_atoms = aseio.read("h2cu_start.in")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)
ase_atoms_jl = ase_atoms.copy()
ase_atoms_jl.pop(-1)
fn = Random.randstring(40) * ".xyz"
aseio.write(fn, ase_atoms_jl)
atoms_julip = read_extxyz(fn)
rm(fn)

model_ml = ace_model("ace_dens_model/h2cu_ace.json", atoms_julip)
density_model = AceLDFA(model_ml; density_unit=u"Å^-3")
model = LDFAFriction(density_model, atoms; friction_atoms=[55, 56])

@testset "ACE LDFA" begin
    F = zeros(3*56, 3*56)
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
