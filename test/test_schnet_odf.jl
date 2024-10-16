ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"  # optional

using Test
using FrictionProviders
using PythonCall
using PyCall: pyimport
using NQCBase: NQCBase
using NQCModels: FrictionModels
using Unitful: @u_str

# SchNet MODEL
spk_utils = pyimport("schnetpack.utils")
spk_envir = pyimport("schnetpack.environment")
spk_ft = pyimport("friction_tensor")
torch = pyimport("torch")
aseio = pyimport("ase.io")

atoms_path = "h2cu_start.in"
model_path = "schnet_model/"
device = "cpu" # vs 'cuda'
cutoff= 5.0

ase_atoms = aseio.read(atoms_path)
atoms, R, cell =  NQCBase.convert_from_ase_atoms(ase_atoms)

model_args = spk_utils.read_from_json("$(model_path)args.json")
model = torch.load("$(model_path)best_model", map_location=device)
environment_provider = spk_envir.AseEnvironmentProvider(cutoff) #loading environment provider, ase in our case, as saved in the args.json file

# setting a calculator 
calculator = spk_ft.FrictionCalculator(model=model, device=device, cutoff=cutoff, friction_tensor="friction_tensor", friction_indices=torch.Tensor([[54,55]]), environment_provider=environment_provider)


schnet_model = SchNetODF(calculator, ase_atoms; friction_unit=u"ps^-1")
odf_model = ODFriction(schnet_model; friction_atoms=[55, 56])

@testset "SchNetModel!" begin
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