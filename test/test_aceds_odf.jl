using Test
using FrictionProviders
using PyCall: pyimport
using NQCBase: NQCBase
using NQCModels: FrictionModels
using Unitful: @u_str
using ACE
using ACEds: ac_matrixmodel
using ACEds.FrictionModels
using ACEds.FrictionModels: Gamma
using JuLIP
using ASE
using JLD2

aseio = pyimport("ase.io")

## LOAD INITIAL ATOMS
atoms_path = "h2cu_start.in"
model_path = "aceds_model/"

ase_atoms = aseio.read(atoms_path)
ase_jl = ASE.ASEAtoms(ase_atoms)
julip_atoms = JuLIP.Atoms(ase_jl)
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)


### LOAD ACEds MODEL

model_params = load_object("$(model_path)h2cu_eft_ace_ds_model_params_m.jld2")

m_cov = ac_matrixmodel(ACE.EuclideanVector(Float64),model_params["species_friction"],model_params["species_env"];n_rep=model_params["n_rep_cov"], rcut_on = model_params["cutoff"], rcut_off = model_params["cutoff"], maxorder_on=model_params["max_order_on"], maxdeg_on=model_params["max_degree_on"],
        species_maxorder_dict_on = model_params["species_maxorder_dict_on_cov"], 
        species_weight_cat_on = model_params["species_weight_cat_on_cov"],
        species_maxorder_dict_off = model_params["species_maxorder_dict_off_cov"], 
        species_weight_cat_off = model_params["species_weight_cat_off_cov"],
        bond_weight = model_params["bond_weight"]
    )
m_equ = ac_matrixmodel(ACE.EuclideanMatrix(Float64),model_params["species_friction"],model_params["species_env"];n_rep=model_params["n_rep_equ"], rcut_on = model_params["cutoff"], rcut_off = model_params["cutoff"], maxorder_on=model_params["max_order_on"], maxdeg_on=model_params["max_degree_on"],
        species_maxorder_dict_on = model_params["species_maxorder_dict_on_equ"], 
        species_weight_cat_on = model_params["species_weight_cat_on_equ"],
        species_maxorder_dict_off = model_params["species_maxorder_dict_off_equ"], 
        species_weight_cat_off = model_params["species_weight_cat_off_equ"],
        bond_weight = model_params["bond_weight"]
    )

c_fit = load_object("$(model_path)h2cu_eft_ace_ds_model_c_fit.jld2")
eft_model_ace = FrictionModel((m_cov,m_equ))
ACE.set_params!(eft_model_ace, c_fit)


### SET THE ODF MODEL
friction_ids = [55,56]
aceds_model = ACEdsODF(eft_model_ace, Gamma, julip_atoms; friction_unit=u"ps^-1")
odf_model = ODFriction(aceds_model; friction_atoms=friction_ids)


### TEST THE MODEL

@testset "ACEdsEFTModel!" begin
    F = zeros(3*56, 3*56)
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(odf_model, F, R)
    F = F[(friction_ids[1]-1)*3+1:friction_ids[2]*3,(friction_ids[1]-1)*3+1:friction_ids[2]*3]

    DoFs = size(R, 1)
    mass_weights = zeros(length(friction_atoms)*DoFs,length(friction_atoms)*DoFs)
    for fx in 1:size(mass_weights,1)
        for fy in 1:size(mass_weights,2)
            mass_weights[fx,fy] = sqrt(atoms.masses[friction_atoms[Int(ceil(fx/DoFs,digits=0))]])*sqrt(atoms.masses[friction_atoms[Int(ceil(fy/DoFs,digits=0))]])
        end
    end

    F ./= mass_weights
    F .= auconvert.(u"ps^-1", F)/u"ps^-1" # F / ps^-1

    F_direct = zeros(length(friction_atoms)*DoFs,length(friction_atoms)*DoFs)
    F_direct .= reinterpret(Matrix,Matrix(Gamma(eft_model_ace, julip_atoms)[friction_atoms, friction_atoms]))

    F ≈ F_direct
end