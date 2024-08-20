using Test
using FrictionProviders
using PyCall: pyimport
using NQCBase: NQCBase
using NQCModels: FrictionModels
using Pandas: read_pickle
using Unitful: @u_str
using ASE
using JuLIP
using ACE1

function ace_model(model_path, cur_atoms)
    IP = ACE1.read_dict(load_dict(model_path)["IP"])
    JuLIP.set_calculator!(cur_atoms, IP)
    
    model = AdiabaticModels.JuLIPModel(cur_atoms)
end

# ACE MODEL
aseio = pyimport("ase.io")

ase_atoms = aseio.read("h2cu_start.in")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)
ase_atoms_jl = ase_atoms.copy()
ase_atoms_jl.pop(-1)
ase_jl = ASE.ASEAtoms(ase_atoms_jl)
atoms_julip = JuLIP.Atoms(ase_jl)

model_ml = ace_model("ace_dens_model/h2cu_ace.json", atoms_julip)
density_model = AceLDFA(model_ml; density_unit=u"Å^-3")
model = LDFAFriction(density_model, atoms; friction_atoms=[55, 56])

@testset "ScikitModel!" begin
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