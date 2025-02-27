using Test
using FrictionProviders
using PythonCall
using NQCBase
using NQCModels: FrictionModels
using Unitful: @u_str

# SCI-KIT MODEL
dscr_d = pyimport("dscribe.descriptors")
aseio = pyimport("ase.io")
pandas = pyimport("pandas")

model_ml = pandas.read_pickle("scikit_models/density_soap_h2cu_lghtst.pkl")
scaler_ml = pandas.read_pickle("scikit_models/scaler_density_soap_h2cu_lghtst.pkl")
ase_atoms = aseio.read("h2cu_start.xyz")
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

desc = dscr_d.SOAP(
    species = PyList(["Cu", "H"]),
    periodic = true,
    r_cut = 3.0, # 7
    n_max = 2, # 12
    l_max = 2, # 8
    average="off" 
)

density_model = SciKitLDFA(desc, model_ml, ase_atoms; density_unit=u"Å^-3", scaler=scaler_ml)
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
