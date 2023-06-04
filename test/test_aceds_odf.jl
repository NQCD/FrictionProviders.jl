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

aceds_model = ACEdsODF(eft_model_ace, Gamma, julip_atoms, JuLIP.set_positions!; friction_unit=u"ps^-1")
odf_model = ODFriction(aceds_model; friction_atoms=[55, 56])


### TEST THE MODEL

@testset "ACEdsModel!" begin
    F = zeros(3*56, 3*56)
    r = @view R[:,1]
    r .= 0
    FrictionModels.friction!(odf_model, F, R)

    for _=1:10
        r .= 0
        r += rand() * cell.vectors[:,1]
        r += rand() * cell.vectors[:,2]
        r += rand() * cell.vectors[:,3]

        FrictionModels.friction!(odf_model, F, R)
    end

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3]
    FrictionModels.friction!(odf_model, F, R)

    r = cell.vectors[:,1] + cell.vectors[:,2] + cell.vectors[:,3] + rand(3)
    FrictionModels.friction!(odf_model, F, R)
end