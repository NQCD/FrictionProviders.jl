"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct ACEdsODF{M,G,A,U}
    "ACE/friction_tensor calculator"
    model::M
    "Gamma ACEds function"
    gamma::G
    "JuLIP atoms object"
    atoms_julip::A
    "Units"
    friction_unit::U
end

function ACEdsODF(model, gamma, atoms_julip; friction_unit=u"ps^-1")
    ACEdsODF(model, gamma, atoms_julip, friction_unit)
end

function set_coordinates!(model::ACEdsODF, R)
    model.atoms_julip.X = ustrip.(auconvert.(u"Å", R'))
end

function friction!(model::ACEdsODF, R::AbstractMatrix, friction::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)
    set_coordinates!(model, R)

    friction = reinterpret(Matrix,Matrix(model.gamma(model.model, model.atoms_julip)[friction_atoms, friction_atoms]))
    friction = austrip.(friction .* model.friction_unit)
end
