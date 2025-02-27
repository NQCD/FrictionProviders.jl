using Test
using FrictionProviders
using PythonCall
using NQCBase: NQCBase
using NQCModels: FrictionModels
using Unitful: @u_str
using UnitfulAtomic
using ACE
using ACEds.FrictionModels
using ACEds.FrictionModels: Gamma
using JuLIP

aseio = pyimport("ase.io")

## LOAD INITIAL ATOMS
atoms_path = "hpt_start.xyz"
model_path = "aceds_model/"

ase_atoms = aseio.read(atoms_path)
julip_atoms = JuLIP.read_extxyz(atoms_path)
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)


### LOAD ACEds MODEL
eft_model_ace = read_dict(load_dict("$(model_path)eft_ac.model"))

### SET THE ODF MODEL
friction_ids = [1]
aceds_model = ACEdsODF(eft_model_ace, Gamma, julip_atoms; friction_unit=u"ps^-1")
odf_model = ODFriction(aceds_model; friction_atoms=friction_ids)


### TEST THE MODEL

@testset "ACEdsEFTModel!" begin
    DoFs = size(R, 1)
    F = zeros(DoFs*length(atoms), DoFs*length(atoms))
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(odf_model, F, R)
    F = F[(friction_ids[1]-1)*DoFs+1:friction_ids[end]*DoFs,(friction_ids[1]-1)*DoFs+1:friction_ids[end]*DoFs]

    mass_weights = zeros(length(friction_ids)*DoFs,length(friction_ids)*DoFs)
    for fx in 1:size(mass_weights,1)
        for fy in 1:size(mass_weights,2)
            mass_weights[fx,fy] = sqrt(atoms.masses[friction_ids[Int(ceil(fx/DoFs,digits=0))]])*sqrt(atoms.masses[friction_ids[Int(ceil(fy/DoFs,digits=0))]])
        end
    end

    F ./= mass_weights
    F .= auconvert.(u"ps^-1", F)/u"ps^-1" # F / ps^-1

    F_direct = zeros(length(friction_ids)*DoFs,length(friction_ids)*DoFs)
    F_direct .= reinterpret(Matrix,Matrix(Gamma(eft_model_ace, julip_atoms)[friction_ids, friction_ids]))

    @test F ≈ F_direct
end
