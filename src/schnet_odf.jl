"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct SchNetODF{C,A,S,U}
    "SchNet/friction_tensor calculator"
    calculator::C
    "Atoms"
    atoms::A
    "Units"
    friction_unit::U
end

function SchNetODF(calculator, atoms_ase; friction_unit=u"ps^-1")
    SchNetODF(calculator, atoms_ase, friction_unit)
end

function set_coordinates!(model::SchNetODF, R)
    model.atoms_ase.set_positions(ustrip.(auconvert.(u"Å", R')))
end

friction!(model.calculator, R, model.friction, model.friction_atoms)

function friction!(model::SchNetODF, R::AbstractMatrix, friction::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Int)
    set_coordinates!(model, R)
    model.calculator.calculate(model.atoms_ase)
    friction = model.calculator.get_friction_tensor()
    friction = austrip(friction .* model.friction_unit)
end
