"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct SchNetODF{C,A,U} <: FrictionModels.TensorialFriction
    "SchNet/friction_tensor calculator"
    calculator::C
    "ASE atoms object"
    atoms_ase::A
    "Units"
    friction_unit::U
    friction_atoms::Union{Vector{Int}, Colon}
    ndofs::Int
end

function SchNetODF(calculator, atoms_ase; friction_unit=u"ps^-1", friction_atoms = 1:size(atoms_ase.positions, 1))
    SchNetODF(calculator, atoms_ase, friction_unit, friction_atoms, size(atoms_ase.positions, 2))
end

function NQCModels.FrictionModels.get_friction_matrix(model::SchNetODF, positions::AbstractMatrix)
    model.atoms_ase.set_positions(au_to_ang.(positions'))
    model.calculator.calculate(model.atoms_ase)
    friction .= model.calculator.get_friction_tensor()
    masses = model.atoms_ase.get_masses()[model.friction_atoms]
    mass_weights = repeat(sqrt.(masses), inner=model.ndofs)' .* repeat(sqrt.(masses), inner=model.ndofs)
    friction .= austrip.(friction .* model.friction_unit)
    friction .*= mass_weights
    friction .= austrip.(friction .* u"u")
end
