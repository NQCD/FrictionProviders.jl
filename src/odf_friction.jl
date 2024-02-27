struct ODFriction{C,D} <: NQCModels.FrictionModels.ElectronicFrictionProvider
    "Friction model"
    eft_model::C
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
    "Cutoff distance"
    cutoff::D
end

ODFriction(eft_model; friction_atoms=collect(range(atoms)), cutoff=Float64(5.0))=ODFriction(eft_model, friction_atoms, cutoff)

NQCModels.ndofs(model::ODFriction) = 3

function FrictionModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)
    indices=friction_matrix_indices(model.friction_atoms, NQCModels.ndofs(model))  
    F[indices, indices] .= get_friction_matrix(model.eft_model, R, model.friction_atoms, model.cutoff) 
end

# Overload friction function to map friction_atoms to Subsystem indices
function FrictionModels.friction!(system::Subsystem{<:ODFriction}, F::AbstractMatrix, R::AbstractMatrix)
    indices=friction_matrix_indices(collect(1:length(system.model.friction_atoms)), NQCModels.ndofs(system))
    F[indices, indices] .= get_friction_matrix(system.model.eft_model, R, system.model.friction_atoms, system.model.cutoff)
end