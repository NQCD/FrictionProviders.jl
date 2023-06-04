"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct SchNetODF{C,A,U}
    "SchNet/friction_tensor calculator"
    calculator::C
    "ASE atoms object"
    atoms_ase::A
    "Units"
    friction_unit::U
end

function SchNetODF(calculator, atoms_ase; friction_unit=u"ps^-1")
    SchNetODF(calculator, atoms_ase, friction_unit)
end

function set_coordinates!(model::SchNetODF, R)
    model.atoms_ase.set_positions(au_to_ang.(R'))
end


function friction!(model::SchNetODF, R::AbstractMatrix, friction::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)
    set_coordinates!(model, R)
    model.calculator.calculate(model.atoms_ase)
    friction .= model.calculator.get_friction_tensor()
    friction = austrip.(friction .* model.friction_unit)
end
