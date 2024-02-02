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

function get_friction_matrix(model::ACEdsODF, R::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)
    set_positions!(model.atoms_julip, au_to_ang.(R))
    DoFs = size(R, 1)
    friction=zeros(eltype(R), length(friction_atoms)*DoFs, length(friction_atoms)*DoFs)
    friction .= reinterpret(Matrix,Matrix(model.gamma(model.model, model.atoms_julip)[friction_atoms, friction_atoms]))
    sqrtmass = sqrt.(model.atoms_julip.M[friction_atoms])
    mass_weights = repeat(sqrtmass * sqrtmass', inner=(DoFs,DoFs))
    @. friction = austrip(friction * model.friction_unit)
    @. friction *= mass_weights
    @. friction = austrip(friction * u"u")
    return friction
end

export get_friction_matrix