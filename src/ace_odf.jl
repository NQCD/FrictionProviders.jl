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

"""
    ACEdsODF(model, gamma, atoms_julip; friction_unit=u"ps^-1")

ACEfriction.jl tensorial electronic friction provider. 

# Arguments

## model
This is the `ACEfriction.FrictionModel` used to predict the friction tensor. 

## gamma
The function used to evaluate the friction tensor from the model using the ACEfriction model and JuLIP Atoms as an input. 

Set this to `ACEfriction.Gamma` unless you need a custom function to postprocess your friction tensor before it is used in dynamics. 

## atoms_julip

A copy of the structure which the friction tensor will be predicted from as `JuLIP.Atoms`. 
This should contain the same atoms and positions, in the same order as your Dynamics Simulation if you are using the model for dynamics.
If you are using isotopes, their atomic masses need to be correctly set within this structure to receive the correct friction. 

## friction_unit

The unit in which the ACEfriction model predicts the friction tensor. By default, this is a relaxation rate tensor in inverse ps. 
 
"""
function ACEdsODF(model, gamma, atoms_julip; friction_unit=u"ps^-1")
    ACEdsODF(model, gamma, atoms_julip, friction_unit)
end

"""
    get_friction_matrix(model::ACEdsODF, R::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)

get_friction_matrix uses an `ACEdsODF` model to predict the friction matrix for `friction_atoms` and return it for just those atoms. 

This behaviour is different to NQCModels.friction!, which returns friction for the whole system, not just `friction_atoms`. 
"""
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
